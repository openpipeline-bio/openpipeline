methods = [
  densmap,
  diffusion_map,
  ivis,
  lmds,
  neuralee,
  pca,
  phate,
  pymde,
  simlr,
  tsne,
  umap
]

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | openproblems_dr_h5mu_to_h5ad.run(
      fromState: [
        "input",
        "input_modality",
        "input_layer_counts",
        "input_layer_normalized",
        "input_var_hvg_score"
      ],
      toState: [
        "method_input": "output"
      ]
    )

    | runEach(
      components: methods,
      filter: { id, state, comp ->
        state.method_id == comp.config.name
      },
      fromState: [
        "input": "method_input"
      ],
      toState: [
        "method_output": "output"
      ]
    )

    | openproblems_dr_h5ad_to_h5mu.run(
      fromState: { id, state ->
        def output_obsm_key = state.output_obsm_key
        if (!output_obsm_key) {
          output_obsm_key = "X_" + state.method_id
        }
        [
          "input_dataset": state.input,
          "input_output": state.method_output,
          "input_modality": state.input_modality,
          "output_obsm_key": output_obsm_key
        ]
      },
      toState: [
        "output": "output"
      ]
    )

    | setState(["output"])

  emit:
  output_ch
}
