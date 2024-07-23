workflow run_wf {
  take:
    input_ch

  main:

    // add id as _meta join id to be able to merge with source channel and end of workflow
    output_ch = input_ch
        // Set aside the output for this workflow to avoid conflicts
        | map {id, state -> 
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
        }
        // Add join_id to state
        | map{ id, state -> 
        def new_state = state + ["_meta": ["join_id": id]]
        [id, new_state]
        }
        | view {"After adding join_id: $it"}

        // Add 'query' id to .obs columns of query dataset
        | add_id.run(
            fromState: {id, state ->
                [
                "input": state.input,
                "input_id": "query",
                "obs_output": "dataset",
                ]
            },
            toState: ["input": "output"])
        // Add 'reference'id to .obs columns of reference dataset
        | add_id.run(
                fromState: {id, state ->
                    [
                    "input": state.reference,
                    "input_id": "reference",
                    "obs_output": "dataset",
                    ]
                },
                toState: ["reference": "output"])
        // Concatenate query and reference datasets
        | concatenate_h5mu.run(
            fromState: { id, state ->
            [
                "input": [state.input, state.reference],
                "input_id": ["query", "reference"],
                "other_axis_mode": "move"
            ]
            },
            toState: ["input": "output"]
            )
        | view {"After concatenation: $it"}
        // Run harmony integration with leiden clustering
        | scgpt_leiden_workflow.run(
            fromState: { id, state ->
            [
                "id": id, //
                "input": state.input,
                "modality": "rna",
                "input_layer": state.input_layer,
                "var_gene_names": state.var_gene_names,
                "obs_batch_label": state.obs_batch_label,
                "model": state.model,
                "model_vocab": state.model_vocab,
                "model_config": state.model_config,
                "finetuned_checkpoints_key": state.finetuned_checkpoints_key,
                "pad_token": state.pad_token,
                "pad_value": state.pad_value,
                "n_hvg": state.n_hvg,
                "max_seq_len": state.max_seq_len,
                "obsm_integrated": state.obsm_integrated,
                "dsbn": state.dsbn,
                "batch_size": state.batch_size,
                "n_input_bins": state.n_input_bins,
                "seed": state.seed, 
                "leiden_resolution": state.leiden_resolution //
            ]
            },
            toState: ["input": "output"]
            )
        | view {"After integration: $it"}
        // Split integrated dataset back into a separate reference and query dataset
        | split_samples.run(
            fromState: { id, state ->
            [
                "input": state.input,
                "modality": "rna",
                "obs_feature": "dataset",
                "output_files": "sample_files.csv",
                "drop_obs_nan": "true",
                "output": "ref_query"
            ]
            },
            toState: [ "output": "output", "output_files": "output_files" ],
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
        // Perform KNN label transfer from reference to query
        | pynndescent_knn.run(
            fromState: { id, state ->
            [
                "input": state.integrated_query,
                "modality": "rna",
                "input_obsm_features": state.obsm_integrated,
                "reference": state.integrated_reference,
                "reference_obsm_features": state.obsm_integrated,
                "reference_obs_targets": state.obs_reference_targets,
                "output_obs_predictions": state.output_obs_predictions,
                "output_obs_probability": state.output_obs_probability,
                "output_compression": state.output_compression,
                "weights": state.weights,
                "n_neighbors": state.n_neighbors,
                "output": state.workflow_output
            ]
            },
            toState: {id, output, state -> ["output": output.output]},
            auto: [ publish: true ]
            )

  emit:
    output_ch
}