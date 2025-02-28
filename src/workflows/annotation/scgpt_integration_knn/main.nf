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
      // Allign query and reference datasets
      | align_query_reference.run(
        fromState: [
          "input": "input",
          "modality": "modality",
          "input_layer": "input_layer",
          "input_var_gene_names": "input_var_gene_names",
          "input_obs_batch": "input_obs_batch_label",
          "reference": "reference",
          "reference_layer": "reference_layer",
          "reference_var_gene_names": "reference_var_gene_names",
          "reference_obs_batch": "reference_obs_batch_label",
          "input_reference_gene_overlap": "input_reference_gene_overlap",
          "overwrite_existing_key": "overwrite_existing_key"
        ],
        args: [
          "input_id": "query",
          "reference_id": "reference",
          "output_layer": "_counts",
          "output_var_gene_names": "_gene_names",
          "output_obs_batch": "_sample_id",
          "output_obs_id": "_dataset"
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
      // Run scgpt integration with leiden clustering
      | scgpt_leiden_workflow.run(
        fromState: { id, state ->
          [
            "id": id,
            "input": state.input,
            "modality": state.modality,
            "model": state.model,
            "model_vocab": state.model_vocab,
            "model_config": state.model_config,
            "finetuned_checkpoints_key": state.finetuned_checkpoints_key,
            "obsm_integrated": state.output_obsm_integrated,
            "pad_token": state.pad_token,
            "pad_value": state.pad_value,
            "n_hvg": state.n_hvg,
            "max_seq_len": state.max_seq_len,
            "dsbn": state.dsbn,
            "batch_size": state.batch_size,
            "n_input_bins": state.n_input_bins,
            "seed": state.seed,
            "leiden_resolution": state.leiden_resolution,
          ]
        },
        args: [
          "input_layer": "_counts",
          "obs_batch_label": "_sample_id",
          "var_gene_names": "_gene_names"
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
           "reference_obs_targets": "reference_obs_targets",
           "output_obs_predictions": "output_obs_predictions",
           "output_obs_probability": "output_obs_probability",
           "output_compression": "output_compression",
           "weights": "weights",
           "n_neighbors": "n_neighbors",
           "output": "workflow_output"
        ],
        toState: {id, output, state -> ["output": output.output]},
      )
    
  emit:
    output_ch
}
