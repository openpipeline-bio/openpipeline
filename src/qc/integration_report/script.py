import scanpy as sc
import mudata as mu
import os
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_integrated_leiden_qc.h5mu",
    "modality": "rna",
    "var_gene_names": None,
    "uns_neighbors": "scGPT_integration_neighbors",
    "obsm_umap": "X_scGPT_umap",
    "obs_batch_label": "sample",
    "obs_cell_label": "cell_type",
    "output": "test_report.pdf",
    "output_raw": True,
    "output_umap_label_fig": "test_umap_label.png",
    "output_umap_batch_fig": "test_umap_batch.png",
    "output_md_report": "test_report.qmd",
    "output_metrics": "test_metrics.json"
    }

meta = {
    "resources_dir": "src/qc/integration_report",
    "functionality_name": "integration_report",
    "temp_dir": "tmp"
}

## VIASH END

# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

logger.info("Reading in data")
# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

# Ensure all batch labels are str
adata.obs[par["obs_batch_label"]] = adata.obs[par["obs_batch_label"]].astype(str)

# Generate vizualizations

logger.info("Generating UMAP visualizations")
fig_cell_type_clusters = sc.pl.embedding(
    adata,
    basis=par["obsm_umap"],
    gene_symbols=par["var_gene_names"],
    neighbors_key=par["uns_neighbors"],
    color=par["obs_cell_label"],
    title="UMAP visualization (label)",
    frameon=False,
    return_fig=True,
    show=False,
)

fig_batch_clusters = sc.pl.embedding(
    adata,
    basis=par["obsm_umap"],
    gene_symbols=par["var_gene_names"],
    neighbors_key=par["uns_neighbors"],
    color=par["obs_batch_label"],
    title="UMAP visualization (batch)",
    frameon=False,
    return_fig=True,
    show=False,
)

# Generate metrics and report

# Open template
with open(meta["resources_dir"] + '/report_template.md', 'r') as file:
    template = file.read()

# Context manager to generate temporary directory to save intermediary files
with tempfile.TemporaryDirectory(
        prefix=f"{meta['functionality_name']}-",
        dir=meta["temp_dir"],
        ignore_cleanup_errors=True) as temp_dir:

    temp_dir = Path(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Directories for intermediary files
    report_md = f"{temp_dir}/report.qmd"
    umap_label_fig = f"{temp_dir}/umap_label.png"
    umap_batch_fig = f"{temp_dir}/umap_batch.png"

    # Save figures
    fig_cell_type_clusters.savefig(umap_label_fig, dpi=300, bbox_inches='tight')
    fig_batch_clusters.savefig(umap_batch_fig, dpi=300, bbox_inches='tight')

    # Fetch variables required for report from the adata object
    logger.info("Fetching integration metrics")
    variables = {
        "overall_integration": round(adata.uns["overall_integration_score"], 2),
        "bio_conservation": round(adata.uns["bio_conservation_score"], 2),
        "batch_correction": round(adata.uns["batch_correction_score"], 2),
        "umap_batch": os.path.basename(umap_batch_fig),
        "umap_label": os.path.basename(umap_label_fig)
    }

    # Round all variables, replace nan's with "-"
    for k, v in adata.uns["bio_conservation_metrics"].items():
        # Check if value is NaN and add to variables dict
        variables[k] = round(v, 2) if v == v else "-"

    for k, v in adata.uns["batch_correction_metrics"].items():
        # Check if value is NaN and add to variables dict
         variables[k] = round(v, 2) if v == v else "-"

    logger.info("Generating report")
    # Fill report template with variables
    markdown_content = template.format(**variables)

    # Save report
    with open(report_md, 'w') as f:
        f.write(markdown_content)

    # PDF rendering
    logger.info("Rendering report to pdf")
    subprocess.run(['quarto', 'render', report_md, '--to', 'pdf'])

    # Move report to output
    logger.info("Publishing of pdf report")
    shutil.move(f"{temp_dir}/report.pdf", par["output"])

    # Save raw data
    if not par["output_raw"]:
        logger.info("Skipping publishing of raw output data because --output_raw was not provided.")
    else:
        logger.info("Publishing of raw output data")
        shutil.move(umap_label_fig, par["output_umap_label_fig"])
        shutil.move(umap_batch_fig, par["output_umap_batch_fig"])
        shutil.move(report_md, par["output_md_report"])

        # Remove figures from variable dict before saving as json
        del variables["umap_batch"]
        del variables["umap_label"]

        metrics = {k: (str(v) if v == v else "-") for k, v in variables.items()}

        with open(par["output_metrics"], 'w') as f:
            json.dump(metrics, f)
