import scib
import mudata as mu
import numpy as np
import warnings

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_integrated_leiden.h5mu",
    "modality": "rna",
    "bio_conservation_metrics": ["nmi","ari","asw_label"],
    "batch_correction_metrics": ["asw_batch","pcr","graph_connectivity"],
    "obsm_embeddings": "X_scGPT",
    "obs_batch_label": "sample",
    "obs_cell_label": "cell_type",
    "uns_neighbors": "scGPT_integration_neighbors",
    "obsp_neighbor_connectivities": "scGPT_integration_connectivities",
    "output_compression": None,
    "obs_cluster": "scGPT_integration_leiden_1.0",
    "output": "resources_test/scgpt/test_resources/Kim2020_Lung_integrated_qc.h5mu",
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

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger.info("Reading in data")
# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

# Remove nan values to calculate scores
if adata.obs["cell_type"].isna().any():
    logger.warning("Cell label obs column contains nan values: removing rows with nan value")
    adata_t = adata[~adata.obs["cell_type"].isna()].copy()
else:
    adata_t = adata.copy()

# Ensure all batch labels are str
adata.obs[par["obs_batch_label"]] = adata.obs[par["obs_batch_label"]].astype(str)

# rename neighbors and connectivities names if necessary
if par["uns_neighbors"] != "neighbors":
    adata_t.uns["neighbors"] = adata_t.uns[par["uns_neighbors"]]
if par["obsp_neighbor_connectivities"] != "connectivities":
    adata_t.obsp["connectivities"] = adata_t.obsp[par["obsp_neighbor_connectivities"]]

# Asses metrics to be calculated
assert len(par["bio_conservation_metrics"]) > 0, "Provide at least one metric for bio conservation."
assert len(par["batch_correction_metrics"]) > 0, "Provide at least one metric for bio conservation."

bio_metric_calc = {
    "nmi": False,
    "ari": False,
    "asw_label": False,
    "isolated_label_f1": False,
    "isolated_label_asw": False,
    "clisi_graph": False,
}

batch_metric_calc = {
    "asw_batch": False,
    "pcr": False,
    "graph_connectivity": False,
    "ilisi_graph": False,
    "kbet": False
}

for bio_metric in par["bio_conservation_metrics"]:
    bio_metric_calc[bio_metric] = True
for batch_metric in par["batch_correction_metrics"]:
    batch_metric_calc[batch_metric] = True

# Batch correction metrics
logger.info(f">> Batch Correction Metrics: {par['batch_correction_metrics']}")
# Dictionary with results of calculated metrics
batch_metrics = {}

# batch_metrics = {}
if batch_metric_calc["asw_batch"]:
    logger.info("Calculating ASW score (batches)")
    asw_batch = scib.metrics.silhouette_batch(
        adata_t,
        batch_key=par["obs_batch_label"],
        label_key=par["obs_cell_label"],
        embed=par["obsm_embeddings"],
        metric="euclidean",
        return_all=False,
        verbose=False,
        )
    batch_metrics["asw_batch"] = asw_batch
else:
    logger.info("Skipping ASW score (batches) calculation")
    batch_metrics["asw_batch"] = np.nan
    
if batch_metric_calc["pcr"]:
    logger.info("Calculating PC regression")
    pcr_score = scib.metrics.pcr_comparison(
        adata_t,
        adata_t,
        embed=par["obsm_embeddings"],
        covariate=par["obs_batch_label"],
        verbose=False
        )
    batch_metrics["pcr"] = pcr_score
else:
    logger.info("Skipping PC Regression calculation")
    batch_metrics["pcr"] = np.nan

if batch_metric_calc["graph_connectivity"]:
    logger.info("Calculating graph connectivity")
    graph_conn_score = scib.metrics.graph_connectivity(
        adata_t,
        label_key=par["obs_cell_label"]
        )
    batch_metrics["graph_connectivity"] = graph_conn_score
else:
    logger.info("Skipping graph connectivity calculation")
    batch_metrics["graph_connectivity"] = np.nan

if batch_metric_calc["ilisi_graph"]:
    logger.info("Calculating graph iLISI")
    graph_ilisi = scib.metrics.lisi.ilisi_graph(
        adata_t,
        batch_key=par["obs_batch_label"],
        type_="knn",
        use_rep=par["obsm_embeddings"]
        )
    batch_metrics["ilisi_graph"] = graph_ilisi
else:
    logger.info("Skipping graph iLISI calculation")
    batch_metrics["ilisi_graph"] = np.nan
    
if batch_metric_calc["kbet"]:
    logger.info("Calculating kBET")
    kbet_score = scib.metrics.kBET(
        adata_t,
        batch_key=par["obs_batch_label"],
        label_key=par["obs_cell_label"],
        type_="knn",
        embed=par["obsm_embeddings"],
        verbose=True
        )
    batch_metrics["kbet"] = kbet_score
else:
    logger.info("Skipping kBET calculation")
    batch_metrics["kbet"] = np.nan
    
# Bio conservation metrics
logger.info(f">>Bio Conservation Metrics: {par['bio_conservation_metrics']}")

# Dictionary to store results of bio conservation metrics 
bio_metrics = {}

if bio_metric_calc["nmi"]:
    logger.info("Calculating NMI score")
    nmi_score = scib.metrics.nmi(
            adata_t,
            cluster_key=par["obs_cluster"],
            label_key=par["obs_cell_label"]
        )
    bio_metrics["nmi"] = nmi_score
else:
    logger.info("Skipping NMI calculation")
    bio_metrics["nmi"] = np.nan
    
if bio_metric_calc["ari"]:
    logger.info("Calculating ARI score")
    ari_score = scib.metrics.ari(
        adata_t,
        cluster_key=par["obs_cluster"],
        label_key=par["obs_cell_label"]
        )
    bio_metrics["ari"] = ari_score
else:
    logger.info("Skipping ARI calculation")
    bio_metrics["ari"] = np.nan
    
if bio_metric_calc["asw_label"]:
    logger.info("Calculating ASW score (cell types)")
    asw_label = scib.metrics.silhouette(
        adata_t,
        label_key=par["obs_cell_label"],
        embed=par["obsm_embeddings"]
        )
    bio_metrics["asw_label"] = asw_label
else:
    logger.info("Skipping ASW score (cell types) calculation")
    bio_metrics["asw_label"] = np.nan
    
if bio_metric_calc["isolated_label_f1"]:
    logger.info("Calculating isolated label F1")
    il_score_f1 = scib.metrics.isolated_labels(
        adata_t,
        label_key=par["obs_cell_label"],
        batch_key=par["obs_batch_label"],
        embed=par["obsm_embeddings"],
        cluster=True
        )
    bio_metrics["isolated_label_f1"] = il_score_f1
else:
    logger.info("Skipping isolated label F1 calculation")
    bio_metrics["isolated_label_f1"] = np.nan

if bio_metric_calc["isolated_label_asw"]:
    logger.info("Calculating isolated label silhouette")
    il_score_asw = scib.metrics.isolated_labels(
        adata_t,
        label_key=par["obs_cell_label"],
        batch_key=par["obs_batch_label"],
        embed=par["obsm_embeddings"],
        cluster=False
        )
    bio_metrics["isolated_label_asw"] = il_score_asw
else:
    logger.info("Skipping isolated label ASW calculation")
    bio_metrics["isolated_label_asw"] = np.nan

if bio_metric_calc["clisi_graph"]:
    logger.info("Calculating graph cLISI")
    graph_clisi = scib.metrics.clisi_graph(
        adata_t,
        label_key=par["obs_cell_label"],
        type_="knn",
        use_rep=par["obsm_embeddings"]
    )
    bio_metrics["clisi_graph"] = graph_clisi
else:
    logger.info("Skipping isolated graph cLISI calculation")
    bio_metrics["clisi_graph"] = np.nan
    
logger.info(">> Aggregate Metrics")
logger.info(f"Calculating Bio Conservation Score based on {par['bio_conservation_metrics']} ")
# Calculate the mean of the bio conservation metrics ignoring the nan values
avg_bio = np.nanmean(list(bio_metrics.values()))

logger.info(f"Calculating Batch Correction Score based on {par['batch_correction_metrics']}")
# Calculate the mean of the bio conservation metrics ignoring the nan values
avg_batch = np.nanmean(list(batch_metrics.values()))

logger.info("Calculating Overall Score")
# Weighted mean of the bio conservation and batch correction scores 
# See https://doi.org/10.1038/s41592-021-01336-8
overall_score = (0.4 * avg_batch) + (0.6 * avg_bio)


logger.info("Writing output data")
adata.uns["bio_conservation_metrics"] = bio_metrics
adata.uns["batch_correction_metrics"] = batch_metrics
adata.uns["bio_conservation_score"] = avg_bio
adata.uns["batch_correction_score"] = avg_batch
adata.uns["overall_integration_score"] = overall_score

mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])
