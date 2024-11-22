---
title: "Integration QC report"
format: pdf
tbl-colwidths: [70,20]
---

## Overview of integration metrics and scores

The overall integration score is **{overall_integration}**. The integration score is a 40/60 weighted average of the batch correction score ({batch_correction}, table 1) and the biological conservation score ({bio_conservation}, table 2).


::: {{layout-ncol=2}}

| metric  | value  |
|--------|--------|
| NMI  | {nmi}   |
| ARI   | {ari}   |
| ASW (label) | {asw_label}   |
| Isolated label (f1) | {isolated_label_f1}   |
| Isolated label (ASW) | {isolated_label_asw}   |
| cLISI   | {clisi_graph}   |
| **biological conservation score** | **{bio_conservation}**   |

: Biological conservation {{.striped}}

| metric  | value  |
|--------|--------|
| ASW (batch)  | {asw_batch}   |
| PCR   | {pcr}   |
| Graph connectivity | {graph_connectivity}   |
| iLISI | {ilisi_graph}   |
| kBET | {kbet}   |
| **batch correction score** | **{batch_correction}**   |

: Batch correction {{.striped}}
:::

The biological conservation metrics consist of NMI (normalized mutual information), ARI (adjusted rand index), ASW (average silhouette width) label, isolated label (f1 and ASW) and graph cLISI (cell-type local inverse Simpson's index) metrics. The batch correction metrics consist of ASW, PCR (principal component regression), graph connectivity, iLISI (integration LISI) and kBET (k-nearest-neighbor batch effect test) metrics.
The biological conservation and batch correction scores are the means of their individual metrics.

All score and metrics calculations are based on the scib (single cell integration benchmark) library (https://doi.org/10.1038/s41592-021-01336-8). 

## Visualizations of the integration embeddings

![](umap_batch.png){{fig-align="center" width=70%}}

![](umap_label.png){{fig-align="center" width=70%}}

The integration embedding figures display the UMAP (uniform manifold approximation and projection) visualizations of the integrated data, colored by batch annotation or cell identity.  