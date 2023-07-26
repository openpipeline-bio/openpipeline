from mudata import read_h5mu
import pandas as pd
import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "uns_key": "metrics_cellranger",
    "output": "foo.h5mu",
    "modality": "rna"
}
### VIASH END


logger.info("Read mudata from file")
mdata = read_h5mu(par['input'])
mod_data = mdata.mod[par['modality']]

logger.info("Joining uns to obs")
# get data frame
uns_df = mod_data.uns[par['uns_key']]

# check for overlapping colnames
intersect_keys = uns_df.keys().intersection(mod_data.obs.keys())
obs_drop = mod_data.obs.drop(intersect_keys, axis=1)

# create data frame to join
uns_df_rep = uns_df.loc[uns_df.index.repeat(mod_data.n_obs)]
uns_df_rep.index = mod_data.obs_names

# create new obs
mod_data.obs = pd.concat([obs_drop, uns_df_rep], axis=1)

logger.info("Write output to mudata file")
mdata.write_h5mu(par['output'], compression=par["output_compression"])

        

