import mudata
import muon as mu
import pandas as pd

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "output_marker_peaks": "marker_peaks.tsv",
    "output_compression": "gzip",
    "modality": "atac",
    "groupby": "leiden",
    "method": "t-test",
}
### VIASH END

def main():
    mdata = mudata.read(par["input"].strip())
    mod = mdata.mod[par["modality"]]

    if "peak_annotation" not in mod.uns["atac"]:
        raise ValueError("Peak annotation not found. Please run `muon.atac.tl.add_peak_annotation` first.")

    mu.atac.tl.rank_peaks_groups(mod, par["groupby"], method=par["method"])

    result = mod.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    marker_peaks_df = pd.DataFrame({
        group + "_" + key: result[key][group]
        for group in groups
        for key in ["names", "genes", "pvals_adj"]
    })

    mdata.mod[par["modality"]].uns["rank_genes_groups"] = mod.uns["rank_genes_groups"]

    mdata.write_h5mu(par["output"].strip(), compression=par["output_compression"])

    marker_peaks_df.to_csv(par["output_marker_peaks"].strip(), sep="\t", index=False)

if __name__ == "__main__":
    main()