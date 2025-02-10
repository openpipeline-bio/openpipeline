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
        // Align query and reference datasets
        | align_query_reference.run(
            fromState: [
                "input": "input",
                "modality": "modality",
                "input_layer": "input_layer",
                "input_obs_batch": "input_obs_batch_label",
                "input_var_gene_names": "input_var_gene_names",
                "reference": "reference",
                "reference_layer": "reference_layer",
                "reference_obs_batch": "reference_obs_batch_label",
                "reference_obs_label": "reference_obs_target",
                "reference_var_gene_names": "reference_var_gene_names",
                "input_reference_gene_overlap": "input_reference_gene_overlap",
                "overwrite_existing_key": "overwrite_existing_key"
            ],
            args: [
                "input_id": "query",
                "reference_id": "reference",
                "output_layer": "_counts",
                "output_var_gene_names": "_gene_names",
                "output_obs_batch": "_sample_id",
                "output_obs_label": "_cell_type",
                "output_obs_id": "_dataset",
                "output_var_common_genes": "_common_vars"
            ],
            toState: [
                "input": "output_query",
                "reference": "output_reference"
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
        | highly_variable_features_scanpy.run(
            fromState: [
                "input": "input",
                "modality": "modality",
                "flavor": "hvg_flavor",
                "n_top_features": "hvg_n_top_features",
            ],
            args: [
                "layer": "_counts",
                "var_input": "_common_vars",
                "var_name_filter": "_common_hvg",
                "obs_batch_key": "_sample_id"
            ],
            toState: [
                "input": "output"
            ]
        )
        | pca.run(
            fromState: [
                "input": "input",
                "modality": "modality",
                "overwrite": "overwrite_existing_key",
                "num_compontents": "pca_num_components"
            ],
            args: [
                "layer": "_counts",
                "var_input": "_common_hvg",
                "obsm_output": "X_pca_harmony",
                "varm_output": "pca_loadings_harmony",
                "uns_output": "pca_variance_harmony",
            ],
            toState: [
                "input": "output"
            ]
        )
        // Run harmony integration with leiden clustering
        | harmony_leiden_workflow.run(
            fromState: { id, state -> [
                "id": id,
                "input": state.input,
                "modality": state.modality,
                "obsm_integrated": state.output_obsm_integrated,
                "theta": state.harmony_theta,
                "leiden_resolution": state.leiden_resolution,
            ]},
            args: [
                "embedding": "X_pca_harmony",
                "uns_neighbors": "harmonypy_integration_neighbors",
                "obsp_neighbor_distances": "harmonypy_integration_distances",
                "obsp_neighbor_connectivities": "harmonypy_integration_connectivities",
                "obs_cluster": "harmony_integration_leiden",
                "obsm_umap": "X_leiden_harmony_umap",
                "obs_covariates": "_sample_id"
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
                "obs_feature": "_dataset",
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
                "reference_obs_targets": "reference_obs_target",
                "output_obs_predictions": "output_obs_predictions",
                "output_obs_probability": "output_obs_probability",
                "output_compression": "output_compression",
                "weights": "knn_weights",
                "n_neighbors": "knn_n_neighbors",
                "output": "workflow_output"
            ],
            toState: {id, output, state -> ["output": output.output]},
        )
    
  emit:
    output_ch
}
