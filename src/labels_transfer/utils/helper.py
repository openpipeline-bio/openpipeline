
def _setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger


def check_arguments(par):
    # check output .obs predictions
    if par["output_obs_predictions"] is None:
        par["output_obs_predictions"] = [ t + "_pred" + t for t in par["reference_obs_targets"]]
    assert len(par["output_obs_predictions"]) == len(par["reference_obs_targets"]), "Number of output_obs_predictions must match number of reference_obs_targets"

    # check output .obs uncertainty
    if par["output_obs_uncertainty"] is None:
        par["output_obs_uncertainty"] = [ t + "_uncertainty" for t in par["reference_obs_targets"]]
    assert len(par["output_obs_uncertainty"]) == len(par["reference_obs_targets"]), "Number of output_obs_uncertainty must match number of reference_obs_targets"

    return 

def get_reference_features(adata_reference, par, logger):
    if par["reference_obsm_key"] is None:
        logger.info("Using .X of reference data")
        train_data = adata_reference.X
    else:
        logger.info(f"Using .obsm[{par['reference_obsm_key']}] of reference data")
        train_data = adata_reference.obsm[par["reference_obsm_key"]]

    return train_data

def get_query_features(adata, par, logger):
    if par["input_obsm_features"] is None:
        logger.info("Using .X of query data")
        query_data = adata.X
    else:
        logger.info(f"Using .obsm[{par['input_obsm_features']}] of query data")
        query_data = adata.obsm[par["input_obsm_features"]]

    return query_data