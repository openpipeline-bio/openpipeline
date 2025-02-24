workflow run_wf {
  take:
    input_ch

  main:
    modalities_ch = input_ch
      // Make sure there is not conflict between the output from this workflow
      // And the output from any of the components
      | map {id, state ->
        def new_state = state + ["workflow_output": state.output]
        [id, new_state]
      }
      // If requested to be detected, make sure the mitochondrial and ribosomal genes
      // are added to the input of the qc metrics calculation
      | map {id, state ->
        def var_qc_default = [state.highly_variable_features_var_output]
        if (state.var_name_mitochondrial_genes) {
          var_qc_default.add(state.var_name_mitochondrial_genes)
        }
        if (state.var_name_ribosomal_genes) {
          var_qc_default.add(state.var_name_ribosomal_genes)
        }
        def newState = state + ["var_qc_metrics": var_qc_default.join(",")]
        [id, newState]
      }
      // If requested, add the id of the events (samples) a column in .obs. 
      // Also allows to make .obs_names (the .obs index) unique, by prefixing the values with an unique id per .h5mu file.
      // The latter is usefull to avoid duplicate observations during concatenation.
      | add_id.run(
        filter: {id, state -> state.add_id_to_obs },
        fromState: {id, state -> 
          def newState = [
            "input": state.input,
            "input_id": id,
            "make_observation_keys_unique": state.add_id_make_observation_keys_unique,
            "obs_output": state.add_id_obs_output,
            "add_id_to_obs": state.add_id_to_obs
          ]
          newState
        },
        toState: {id, output, state ->
          def keysToRemove = ["add_id_to_obs", "add_id_obs_output", "add_id_make_observation_keys_unique"]
          def newState = state.findAll{it.key !in keysToRemove}
          newState + ["input": output.output]
        }
      )
      | split_modalities_workflow.run(
        fromState: {id, state ->
          def newState = ["input": state.input, "id": id]
        },
        toState: ["output": "output", "output_types": "output_types"]
      )
      | flatMap {id, state ->
        def outputDir = state.output
        def types = readCsv(state.output_types.toUriString())
        
        types.collect{ dat ->
          // def new_id = id + "_" + dat.name
          def new_id = id // it's okay because the channel will get split up anyways
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["input": new_data, modality: dat.name]]
        }
      }
      // Remove arguments from split modalities from state
      | map {id, state -> 
        def keysToRemove = ["output_types"]
        def newState = state.findAll{it.key !in keysToRemove}
        [id, newState]
      }
      | view {"After splitting modalities: $it"}


    //
    // Singlesample processing
    //
    def singlesample_arguments = [
      "rna": [
        "min_counts": "rna_min_counts",
        "max_counts": "rna_max_counts",
        "min_genes_per_cell": "rna_min_genes_per_cell",
        "max_genes_per_cell": "rna_max_genes_per_cell",
        "min_cells_per_gene": "rna_min_cells_per_gene",
        "min_fraction_mito": "rna_min_fraction_mito",
        "max_fraction_mito": "rna_max_fraction_mito",
        "var_name_mitochondrial_genes": "var_name_mitochondrial_genes",
        "obs_name_mitochondrial_fraction": "obs_name_mitochondrial_fraction",
        "var_gene_names": "var_gene_names",
        "mitochondrial_gene_regex": "mitochondrial_gene_regex",
        "min_fraction_ribo": "rna_min_fraction_ribo",
        "max_fraction_ribo": "rna_max_fraction_ribo",
        "var_name_ribosomal_genes": "var_name_ribosomal_genes",
        "obs_name_ribosomal_fraction": "obs_name_ribosomal_fraction",
        "ribosomal_gene_regex": "ribosomal_gene_regex",
        "layer": "rna_layer"
      ],
      "prot": [
        "min_counts": "prot_min_counts",
        "max_counts": "prot_max_counts",
        "min_proteins_per_cell": "prot_min_proteins_per_cell",
        "max_proteins_per_cell": "prot_max_proteins_per_cell",
        "min_cells_per_protein": "prot_min_cells_per_protein",
        "layer": "prot_layer",
      ],
      "gdo": [
        "min_counts": "gdo_min_counts",
        "max_counts": "gdo_max_counts",
        "min_guides_per_cell": "gdo_min_guides_per_cell",
        "max_guides_per_cell": "gdo_max_guides_per_cell",
        "min_cells_per_guide": "gdo_min_cells_per_guide",
        "layer": "gdo_layer",
      ], 
    ].asImmutable()

    multisample_ch_known = modalities_ch 
      // run the singlesample processing
      | runEach(
        components: [rna_singlesample, prot_singlesample, gdo_singlesample],
        filter: { id, state, component ->
          state.modality + "_singlesample" == component.config.name
        },
        fromState: { id, state, component ->
          def newState = singlesample_arguments.get(state.modality).collectEntries{key_, value_ -> 
            [key_, state[value_]]
          }
          return newState + ["id": id, "input": state.input]
        },
        toState: ["input": "output"],
      )

    multisample_ch_unknown = modalities_ch
      | filter{id, state -> state.modality !in singlesample_arguments.keySet()}

    output_ch = multisample_ch_unknown.mix(multisample_ch_known)
      // Remove arguments for singlesample processing from state.
      | map {id, state -> 
        def keysToRemove = singlesample_arguments.inject([]){currentKeys, modality, stateMapping -> 
            currentKeys += stateMapping.values()
        }
        def allwayskeep = ["gdo_layer", "rna_layer", "prot_layer", "workflow_output"]
        def newState = state.findAll{(it.key !in keysToRemove + ["id"]) || (it.key in allwayskeep)}
        [id, newState]
      }
      | view {"After singlesample processing: $it"}

      //
      // Concatenation: join observations across samples together per modality. 
      //
      // Concatenate multiple single-sample unimodal MuData objects back into several multi-sample files.
      // One multi-sample MuData file is created per modality.
      //
      | map { id, state -> // Put modality name in first element so that we can group on it
        [state.modality, id, state]
      }
      | groupTuple(by: 0, sort: "hash")
      | view {"After groupTuple: $it"}
      | map { modality, old_ids, states ->
        def new_id = "combined_$modality"
        // keys in the new state that should not have a unique value across samples
        def new_state_non_unique_values = [
          "input": states.collect{it.input},
          "input_id": old_ids,
          "_meta": ["join_id": old_ids[0]]
        ]
        // Gather the keys from the different states,
        // one state might contain more keys compared to another (so create a set)
        def all_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
        }.minus(["output", "input_id", "input", "_meta"])
        // Create the new state from the keys, values should be the same across samples
        def new_state = all_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across samples. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            // take the unique value from the set (there is only one)
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
        }
        def final_state = new_state_non_unique_values + new_state
        [new_id, final_state]
      }
      | concatenate_h5mu.run(
        fromState: [
          "input": "input",
          "input_id": "input_id"
        ],
        toState: {id, output, state -> 
          def keysToRemove = ["input_id"]
          def newState = state.findAll{it.key !in keysToRemove}
          newState + ["input": output.output]
        }, 
      )

      | view {"After concatenation: $it"}
      | toSortedList()
      | map {modalities_states -> 
        def states = modalities_states.collect{it[1]}
        def new_input = states.collect{it.input}
        def join_id = states[0]._meta.join_id
        def other_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
          }.minus(["output", "input", "modality", "_meta"])
        def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across modalities. Please report this \
                                                 as a bug. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
          }
         ["merged", new_state + ["input": new_input, "_meta": ["join_id": join_id]]]
      }
      | process_batches.run(
        fromState: {id, state ->
          [
            "id": id,
            "input": state.input,
            "output": state.workflow_output,
            "highly_variable_features_var_output": state.highly_variable_features_var_output,
            "highly_variable_features_obs_batch_key": state.highly_variable_features_obs_batch_key,
            "var_qc_metrics": state.var_qc_metrics,
            "top_n_vars": state.top_n_vars, 
            "pca_overwrite": state.pca_overwrite,
            "rna_layer": state.rna_layer,
            "prot_layer": state.prot_layer,
            "clr_axis": state.clr_axis,
            "rna_enable_scaling": state.rna_enable_scaling,
            "rna_scaling_output_layer": state.rna_scaling_output_layer,
            "rna_scaling_pca_obsm_output": state.rna_scaling_pca_obsm_output,
            "rna_scaling_pca_loadings_varm_output": state.rna_scaling_pca_loadings_varm_output,
            "rna_scaling_pca_variance_uns_output": state.rna_scaling_pca_variance_uns_output,
            "rna_scaling_umap_obsm_output": state.rna_scaling_umap_obsm_output,
            "rna_scaling_max_value": state.rna_scaling_max_value,
            "rna_scaling_zero_center": state.rna_scaling_zero_center,
          ]
        },
        toState: {id, output, state -> 
          [
            "output": output.output,
            "_meta": state._meta,
          ]
        }
      )
      | view {"After process_batches: $it"}

  emit:
    output_ch
}
