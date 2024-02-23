import mudata as mu
from pathlib import Path
from operator import attrgetter, itemgetter
from pandas import Series
import re
import numpy as np


### VIASH START
par = {
    "input": "./resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
    "modality": "rna",
    "matrix": "var",
    "input_column": "gene_symbol",
    "regex_pattern": "^[mM][tT]-",
    "output": "foo.h5mu",
    "input_id": "mouse",
    "output_match_column": "test",
    "output_fraction_column": "fraction_test",
    "output_compression": "gzip"
}
### VIASH END

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

def main(par):
    input_file, output_file, mod_name = Path(par["input"]), Path(par["output"]), par['modality']
    logger.info(f"Compiling regular expression '{par['regex_pattern']}'.")
    try:
        compiled_regex = re.compile(par["regex_pattern"])
    except (TypeError, re.error) as e:
        raise ValueError(f"{par['regex_pattern']} is not a valid regular expression pattern.") from e
    else:
        if compiled_regex.groups:
            raise NotImplementedError("Using match groups is not supported by this component.")
    logger.info('Reading input file %s, modality %s.', input_file, mod_name)

    mudata = mu.read_h5mu(input_file)
    modality_data = mudata[mod_name]
    logger.info("Reading input file done.")
    logger.info("Using annotation dataframe '%s'.", par["matrix"])
    annotation_matrix = getattr(modality_data, par['matrix'])
    default_column = {
        "var": attrgetter("var_names"),
        "obs": attrgetter("obs_names")
    }
    if par["input_column"]:
        logger.info("Input column '%s' was specified.", par["input_column"])
        try:
            annotation_column = annotation_matrix[par["input_column"]]
        except KeyError as e:
            raise ValueError(f"Column {par['input_column']} could not be found for modality "
                            f"{par['modality']}. Available columns:"
                            f" {','.join(annotation_matrix.columns.to_list())}") from e
    else:
        logger.info(f"No input column specified, using '.{par['matrix']}_names'")
        annotation_column = default_column[par['matrix']](modality_data).to_series()
    logger.info("Applying regex search.")
    grep_result = annotation_column.str.contains(par["regex_pattern"], regex=True)
    logger.info("Search results: %s", grep_result.value_counts())

    other_axis_attribute = {
        "var": "obs",
        "obs": "var"
    }
    if par['output_fraction_column']:
        logger.info("Enabled writing the fraction of values that matches to the pattern.")
        input_layer = modality_data.X if not par["input_layer"] else modality_data.layers[par["input_layer"]]
        pct_matching = np.ravel(np.sum(input_layer[:, grep_result], axis=1) / np.sum(input_layer, axis=1))
        assert ((pct_matching >= 0) & (pct_matching <= 1)).all(), \
                "Fractions are not within bounds, please report this as a bug"
        logger.info("Fraction statistics: \n%s", Series(pct_matching).describe())
        output_matrix = other_axis_attribute[par['matrix']]
        logger.info("Writing fractions to matrix '%s', column '%s'",
                    output_matrix, par['output_fraction_column'])
        getattr(modality_data, output_matrix)[par['output_fraction_column']] = pct_matching
    logger.info("Adding values that matched the pattern to '%s', column '%s'",
                par["matrix"], par["output_match_column"])
    getattr(modality_data, par['matrix'])[par["output_match_column"]] = grep_result
    logger.info("Writing out data to '%s' with compression '%s'.",
                output_file, par["output_compression"])
    mudata.write(output_file, compression=par["output_compression"])

if __name__ == "__main__":
    main(par)