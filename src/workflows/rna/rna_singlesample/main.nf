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
    // Check for correctness of mitochondrial and ribosomal gene detection arguments
    | map { id, state ->
      def new_state = [:]
      if (state.obs_name_mitochondrial_fraction && !state.var_name_mitochondrial_genes) {
        throw new RuntimeException("Using --obs_name_mitochondrial_fraction requires --var_name_mitochondrial_genes.")
      }
      if (!state.obs_name_mitochondrial_fraction && state.var_name_mitochondrial_genes) {
        new_state.obs_name_mitochondrial_fraction = "fraction_${state.var_name_mitochondrial_genes}"
      }
      if ((state.min_fraction_mito != null || state.max_fraction_mito != null) && !state.var_name_mitochondrial_genes) {
        throw new RuntimeException("Enabling --min_fraction_mito or --max_fraction_mito requires --var_name_mitochondrial_genes.")
      }
      if (state.var_gene_names && !state.var_name_mitochondrial_genes && !state.var_name_ribosomal_genes) {
        System.err.println("Warning: --var_gene_names requires either --var_name_mitochondrial_genes or --var_name_ribosomal_genes to be set.")
      }
      if (state.mitochondrial_gene_regex && !state.var_name_mitochondrial_genes && !state.var_name_ribosomal_genes) {
        System.err.println("Warning: --mitochondrial_gene_regex requires either --var_name_mitochondrial_genes or --var_name_ribosomal_genes to be set.")
      }
      if (state.obs_name_ribosomal_fraction && !state.var_name_ribosomal_genes) {
        throw new RuntimeException("Using --obs_name_ribosomal_fraction requires --var_name_ribosomal_genes.")
      }
      if (!state.obs_name_ribosomal_fraction && state.var_name_ribosomal_genes) {
        new_state.obs_name_ribosomal_fraction = "fraction_${state.var_name_ribosomal_genes}"
      }
      if ((state.min_fraction_ribo != null || state.max_fraction_ribo != null) && !state.var_name_ribosomal_genes) {
        throw new RuntimeException("Enabling --min_fraction_ribo or --max_fraction_ribo requires --var_name_ribosomal_genes.")
      }
      [id, state + new_state]
    }
    | qc.run(
      fromState: { id, state ->
        // The rna singlesample processing allows detecting mitochondrial genes and filtering based
        // on the fraction of mitochondrial genes per cell
        // This behaviour is optional based on the presence of var_name_mitochondrial_genes
        // The behavior of other components must be tuned to this argument as well
        def args = [
          "id": id,
          "input": state.input,
          // disable other qc metric calculations
          // only mitochondrial gene detection is required at this point
          "top_n_vars": [],
          "output_obs_num_nonzero_vars": null,
          "output_obs_total_counts_vars": null,
          "output_var_num_nonzero_obs": null,
          "output_var_total_counts_obs": null,
          "output_var_obs_mean": null,
          "output_var_pct_dropout": null,
          "output": state.output,
          "modality": "rna",
          "layer": state.layer,
        ]
        
        if (state.var_name_mitochondrial_genes) {
          // Check if user has defined var columns to calculate metrics
          def new_var_qc_metrics = state.var_qc_metrics != null ? state.var_qc_metrics : []
          assert new_var_qc_metrics instanceof List
          // Add the mitochondrial genes var column to the columns to calculate statistics for if set.
          new_var_qc_metrics = ((new_var_qc_metrics as Set) + [state.var_name_mitochondrial_genes]) as List
          
          args += [
            "var_qc_metrics": new_var_qc_metrics,
            "obs_name_mitochondrial_fraction": state.obs_name_mitochondrial_fraction,
            "var_gene_names": state.var_gene_names,
            "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
            "mitochondrial_gene_regex": state.mitochondrial_gene_regex
          ]
        }

        if (state.var_name_ribosomal_genes) {
          // Check if user has defined var columns to calculate metrics
          def new_var_qc_metrics = state.var_qc_metrics != null ? state.var_qc_metrics : []
          assert new_var_qc_metrics instanceof List
          // Add the ribosomal genes var column to the columns to calculate statistics for if set.
          new_var_qc_metrics = ((new_var_qc_metrics as Set) + [state.var_name_ribosomal_genes]) as List
          
          args += [
            "var_qc_metrics": new_var_qc_metrics,
            "obs_name_ribosomal_fraction": state.obs_name_ribosomal_fraction,
            "var_gene_names": state.var_gene_names,
            "var_name_ribosomal_genes": state.var_name_ribosomal_genes,
            "ribosomal_gene_regex": state.ribosomal_gene_regex
          ]
        }

        return args
      },
      toState: ["input": "output"]
    )
    | delimit_fraction.run(
      runIf: {id, state -> state.var_name_mitochondrial_genes},
      fromState: {id, state -> 
      [
        "input": state.input,
        "obs_name_filter": "filter_mitochondrial",
        "min_fraction": state.min_fraction_mito,
        "max_fraction": state.max_fraction_mito,
        "obs_fraction_column": state.obs_name_mitochondrial_fraction,
      ]
      },
      toState: ["input": "output"]
    )
    | delimit_fraction.run(
      runIf: {id, state -> state.var_name_ribosomal_genes},
      fromState: {id, state -> 
      [
        "input": state.input,
        "obs_name_filter": "filter_ribosomal",
        "min_fraction": state.min_fraction_ribo,
        "max_fraction": state.max_fraction_ribo,
        "obs_fraction_column": state.obs_name_ribosomal_fraction,
      ]
      },
      toState: ["input": "output"]
    )
    // cell filtering
    | filter_with_counts.run(
      key: "rna_filter_with_counts",
      fromState: { id, state ->
        [
          "input": state.input,
          "layer": state.layer,
          "obs_name_filter": "filter_with_counts",
          "var_name_filter": "filter_with_counts",
          "min_counts": state.min_counts,
          "max_counts": state.max_counts,
          "min_genes_per_cell": state.min_genes_per_cell,
          "max_genes_per_cell": state.max_genes_per_cell,
          "min_cells_per_gene": state.min_cells_per_gene,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "rna_do_filter",
      fromState: {id, state ->
        // do_filter does not need a layer argument because it filters all layers
        // from a modality.
        def stateMapping = [
          input: state.input,
          var_filter: ["filter_with_counts"]
        ]
        def obs_filter = ["filter_with_counts"]
        if (state.var_name_mitochondrial_genes) {
          obs_filter += ["filter_mitochondrial"]
        }
        if (state.var_name_ribosomal_genes) {
          obs_filter += ["filter_ribosomal"]
        }
        stateMapping += ["obs_filter": obs_filter]
        return stateMapping
      },
      toState: ["input": "output"]
    )
    // doublet calling
    | filter_with_scrublet.run(
      fromState: [
        "input": "input",
        "output": "workflow_output",
        "layer": "layer",
      ],
      args: [output_compression: "gzip"],
    )
    | setState(["output": "output"])

  emit:
  output_ch
}