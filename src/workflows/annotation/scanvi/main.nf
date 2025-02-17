workflow run_wf {
  take:
    input_ch

  main:


    output_ch = input_ch
        // Set aside the output for this workflow to avoid conflicts
        | map {id, state -> 
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
        }
        // add id as _meta join id to be able to merge with source channel and end of workflow
        | map{ id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        [id, new_state]
        }
        | view {"After adding join_id: $it"}
        // Add 'query' id to .obs columns of query dataset
        | add_id.run(
            fromState: [
                "input": "input",
            ],
            args:[
                "input_id": "query",
                "obs_output": "dataset",
            ],
            toState: ["input": "output"]
        )
        // Add 'reference'id to .obs columns of reference dataset
        | add_id.run(
                fromState:[
                    "input": "reference",
                ],
                args:[
                    "input_id": "reference",
                    "obs_output": "dataset"
                ],
                toState: ["reference": "output"]
        )
        // Make sure that query and reference dataset have batch information in the same .obs column
        // By copying the respective .obs columns to the obs column "batch_label"
        | copy_obs.run(
            fromState: [
                "input": "input",
                "modality": "modality",
                "input_obs_key": "input_obs_batch_label",
            ],
            args: [
                "output_obs_key": "batch_label"
            ],
            toState: [
                "input": "output"
            ]
        )
        | copy_obs.run(
            fromState: [
                "input": "reference",
                "modality": "modality",
                "input_obs_key": "reference_obs_batch_label",
            ],
            args: [
                "output_obs_key": "batch_label"
            ],
            toState: [
                "reference": "output"
            ]
        )
        // Concatenate query and reference datasets prior to integration
        | concatenate_h5mu.run(
            fromState: { id, state -> [
                "input": [state.input, state.reference]
                ]
            },
            args: [
                "input_id": ["query", "reference"],
                "other_axis_mode": "move"
            ],
            toState: ["input": "output"]
        )
        | view {"After concatenation: $it"}
        // Run scvi integration with leiden clustering
        | scvi_leiden_workflow.run(
            fromState: { id, state ->
            [
                "id": id,
                "input": state.input,
                "layer": state.layer, 
                "modality": state.modality,
                "obsm_output": state.output_obsm_integrated,
                "leiden_resolution": state.leiden_resolution,
                "var_input": state.var_hvg,
                "early_stopping": state.early_stopping,
                "early_stopping_monitor": state.early_stopping_monitor,
                "early_stoping_patience": state.early_stoping_patience,
                "early_stopping_min_delta": state.early_stopping_min_delta,
                "max_epochs": state.max_epochs,
                "reduce_lr_on_plateau": state.reduce_lr_on_plateau,
                "lr_factor": state.lr_factor,
                "lr_patience": state.lr_patience
            ]},
            args: [
                "uns_neighbors": "scvi_integration_neighbors",
                "obsp_neighbor_distances": "scvi_integration_distances",
                "obsp_neighbor_connectivities": "scvi_integration_connectivities",
                "obs_cluster": "scvi_integration_leiden",
                "obsm_umap": "X_leiden_scvi_umap",
                "obs_batch": "batch_label"
            ],
            toState: ["input": "output"]
        )
        | view {"After integration: $it"}
        // Split integrated dataset back into a separate reference and query dataset
        | split_h5mu.run(
            fromState: [
                "input": "input",
                "modality": "modality"
            ],
            args: [
                "obs_feature": "dataset",
                "output_files": "sample_files.csv",
                "drop_obs_nan": "true",
                "output": "ref_query"
            ],
            toState: [ 
                "output": "output", 
                "output_files": "output_files" 
            ],
            auto: [ publish: true ]
        )
        | view {"After sample splitting: $it"}
        // map the integrated query and reference datasets back to the state
        | map {id, state ->
            def outputDir = state.output
            def files = readCsv(state.output_files.toUriString())
            def query_file = files.findAll{ dat -> dat.name == 'query' }
            assert query_file.size() == 1, 'there should only be one query file'
            def reference_file = files.findAll{ dat -> dat.name == 'reference' }
            assert reference_file.size() == 1, 'there should only be one reference file'
            def integrated_query = outputDir.resolve(query_file.filename)
            def integrated_reference = outputDir.resolve(reference_file.filename)
            def newKeys = ["integrated_query": integrated_query, "integrated_reference": integrated_reference]
            [id, state + newKeys]
            }
        | view {"After splitting query: $it"}
        // Perform KNN label transfer from integrated reference to integrated query
        | knn.run(
            fromState: [
                "input": "integrated_query",
                "modality": "modality",
                "input_obsm_features": "output_obsm_integrated",
                "reference": "integrated_reference",
                "reference_obsm_features": "output_obsm_integrated",
                "reference_obs_targets": "reference_obs_targets",
                "output_obs_predictions": "output_obs_predictions",
                "output_obs_probability": "output_obs_probability",
                "output_compression": "output_compression",
                "weights": "weights",
                "n_neighbors": "n_neighbors",
                "output": "workflow_output"
            ],
            toState: {id, output, state -> ["output": output.output]},
            auto: [ publish: true ]
        )

  emit:
    output_ch
}