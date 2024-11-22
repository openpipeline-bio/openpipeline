workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Calculate integration qc metrics and add them to the MuData file
    | integration_metrics.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "obsm_output": "obsm_output",
        "bio_conservation_metrics": "bio_conservation_metrics",
        "batch_correction_metrics": "batch_correction_metrics",
        "obsm_embeddings": "obsm_embeddings",
        "obs_batch_label": "obs_batch_label",
        "obs_cell_label": "obs_cell_label",
        "obs_cluster": "obs_cluster",
        "uns_neighbors": "uns_neighbors",
        "obsp_neighbor_connectivities": "obsp_neighbor_connectivities",
      ],
      toState: [
        "input": "output"
      ]
    )
    
    // Generate visualizations and a report based on the metrics
    | integration_report.run(
      fromState: [
        "input": "input",
        "modality": "modality",
        "var_gene_names": "var_gene_names",
        "obs_batch_label": "obs_batch_label",
        "obs_cell_label": "obs_cell_label",
        "uns_neighbors": "uns_neighbors",
        "obsm_umap": "obsm_umap",
        "output_raw": "output_raw"
      ],
      toState: { id, output, state ->
        output
      }
    )

  emit:
  output_ch
}