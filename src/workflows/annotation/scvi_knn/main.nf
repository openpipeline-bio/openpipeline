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
      | view {"After adding join_id: $it"}
      // Allign query and reference datasets
      | align_query_reference.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "input_layer": "input_layer",
          "input_layer_lognormalized": "input_layer_lognormalized",
          "input_obs_batch": "input_obs_batch_label",
          "input_var_gene_names": "input_var_gene_names",
          "reference": "reference",
          "reference_layer": "reference_layer",
          "reference_layer_lognormalized": "reference_layer_lognormalized",
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
          "output_layer_lognormalized": "_log_normalized",
          "output_var_gene_names": "_gene_names",
          "output_obs_batch": "_sample_id",
          "output_obs_label": "_cell_type",
          "output_obs_id": "_dataset",
          "output_var_common_genes": "_common_vars",
          "align_layers_lognormalized_counts": "true"
        ],
        toState: [
          "input": "output_query",
          "reference": "output_reference"
        ]
      )
      // Concatenate query and reference datasets prior to integration
      | concatenate_h5mu.run(
        fromState: { id, state -> 
          ["input": [state.input, state.reference]]
        },
        args: [
          "input_id": ["query", "reference"],
          "other_axis_mode": "move"
        ],
        toState: ["input": "output"]
      )
      | view {"After concatenation: $it"}
      // Calculate HVG across query and reference
      | highly_variable_features_scanpy.run(
          fromState: [
            "input": "input",
            "modality": "modality",
            "n_top_features": "n_hvg"

          ],
          args: [
            "layer": "_log_normalized",
            "var_input": "_common_vars",
            "var_name_filter": "_common_hvg",
            "obs_batch_key": "_sample_id"
          ],
          toState: [
            "input": "output"
          ]
        )
      | delete_layer.run(
          key: "delete_aligned_lognormalized_counts_layer",
          fromState: [
            "input": "input",
            "modality": "modality",
          ],
          args: [
            "layer": "_log_normalized",
            "missing_ok": "true"
          ],
          toState: [
            "input": "output"
          ]
      )
      // Run scvi integration with leiden clustering
      | scvi_leiden_workflow.run(
        fromState: { id, state -> [
          "id": id,
          "input": state.input,
          "modality": state.modality,
          "obsm_output": state.output_obsm_integrated,
          "leiden_resolution": state.leiden_resolution,
          "early_stopping": state.scvi_early_stopping,
          "early_stopping_monitor": state.scvi_early_stopping_monitor,
          "early_stoping_patience": state.scvi_early_stoping_patience,
          "early_stopping_min_delta": state.scvi_early_stopping_min_delta,
          "max_epochs": state.scvi_max_epochs,
          "reduce_lr_on_plateau": state.scvi_reduce_lr_on_plateau,
          "lr_factor": state.scvi_lr_factor,
          "lr_patience": state.scvi_lr_patience
        ]},
        args: [
          "var_input": "_common_hvg",
          "uns_neighbors": "scvi_integration_neighbors",
          "obsp_neighbor_distances": "scvi_integration_distances",
          "obsp_neighbor_connectivities": "scvi_integration_connectivities",
          "obs_cluster": "scvi_integration_leiden",
          "obsm_umap": "X_leiden_scvi_umap",
          "obs_batch": "_sample_id",
          "layer": "_counts"
        ],
        toState: ["input": "output"]
      )
      | delete_layer.run(
          key: "delete_aligned_counts_layer",
          fromState: [
            "input": "input",
            "modality": "modality",
          ],
          args: [
            "layer": "_counts",
            "missing_ok": "true"
          ],
          toState: [
            "input": "output"
          ]
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
           "weights": "knn_weights",
           "n_neighbors": "knn_n_neighbors",
           "output": "workflow_output"
        ],
        args:[
          "output_compression": "gzip"
        ],
        toState: {id, output, state -> ["output": output.output]}
      )
    
  emit:
    output_ch
}
