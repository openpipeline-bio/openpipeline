import scanpy as sc
import muon as mu
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": ["rna"],
    "output_key": "pca",
    "num_components": 25,
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

okey = par["output_key"]

logger.info("Reading %s.", par["input"])
mdata = mu.read_h5mu(par["input"])

for mod in par['modality']:
    logger.info("Computing PCA components for modality '%s'", mod)
    data = mdata.mod[mod]
    if par['layer'] and par['layer'] not in data.layers:
        raise ValueError(f"{par['layer']} was not found in modality {mod}.")
    layer = data.X if not par['layer'] else data.layers[par['layer']]
    # run pca
    # sc.tl.pca(data, n_comps=par["num_components"])
    X_pca, loadings, variance, variance_ratio = sc.tl.pca(
        layer, 
        n_comps=par["num_components"], 
        return_info=True
    )

    layer_name = "X_" if not par['layer'] else f"{par['layer']}_"
    # store output in specific objects
    data.obsm[layer_name+okey] = X_pca
    data.varm["loadings_"+okey] = loadings.T
    data.uns[okey] = { "variance": variance, "variance_ratio": variance_ratio }

logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"])

logger.info("Finished")