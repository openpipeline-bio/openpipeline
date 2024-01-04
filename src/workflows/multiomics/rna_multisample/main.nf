workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    | normalize_total.run(
      fromState: { id, state ->
        [
          "input": state.input,
          "output_layer": "normalized",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"],
    )
    | log1p.run( 
      fromState: { id, state ->
        [
          "input": state.input,
          "output_layer": "log_normalized",
          "input_layer": "normalized",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"]
    )
    | delete_layer.run(
      fromState: {id, state -> 
        [
          "input": state.input,
          "layer": "normalized",
          "modality": state.modality
        ]
      },
      toState: ["input": "output"]
    )
    | filter_with_hvg.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "layer": "log_normalized",
          "modality": state.modality,
          "var_name_filter": state.filter_with_hvg_var_output,
          "n_top_genes": state.filter_with_hvg_n_top_genes,
          "flavor": state.filter_with_hvg_flavor,
          "obs_batch_key": state.filter_with_hvg_obs_batch_key
        ]
      },
      toState: ["input": "output"],
    )
    | rna_qc.run(
      // TODO: remove when viash 0.8.3 is released
      key: "rna_qc",
      // layer: null to use .X and not log transformed
      fromState: {id, state ->
        [
          "id": id,
          "input": state.input,
          "output": state.workflow_output,
          "input_layer": null,
          "output_compression": "gzip",
          "modality": state.modality,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
          "num_nonzero_vars": state.num_nonzero_vars,
          "total_counts_var": state.total_counts_var,
          "num_nonzero_obs": state.num_nonzero_obs,
          "total_counts_obs": state.total_counts_obs,
          "obs_mean": state.obs_mean,
          "pct_dropout": state.pct_dropout
        ]
      },
    )
    | setState(["output"])

  emit:
  output_ch
}