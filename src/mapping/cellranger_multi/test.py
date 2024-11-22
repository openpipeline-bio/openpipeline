from pathlib import Path
import tarfile
from textwrap import dedent
import re
import pytest
import sys
import subprocess

## VIASH START
meta = {
    'executable': './target/docker/mapping/cellranger_multi/cellranger_multi',
    'resources_dir': 'resources_test/',
    'cpus': 15,
    'memory_gb': 20,
    'config': 'src/mapping/cellranger_multi/config.vsh.yaml'
}
## VIASH END


def make_path_relative(some_path):
    absolute_input_path = Path(some_path).resolve()
    absolute_cwd = Path.cwd().resolve()
    try:
        return absolute_input_path.relative_to(absolute_cwd)
    except ValueError as e:
    # TODO: python 3.12: remove lines below and add walk_up=True to `relative_to` call
        if "is not in the subpath of" in str(e):
            _, *parts_input = absolute_input_path.parts
            _, *parts_cwd = absolute_cwd.parts
            parts_input.reverse()
            parts_cwd.reverse()
            while parts_cwd and parts_input and parts_cwd[-1] == parts_input[-1]:
                parts_input.pop()
                parts_cwd.pop()
            for part in parts_cwd:
                if not part or part == '.':
                    pass
                else:
                    parts_input.append('..')
            relative_path = type(absolute_input_path)('', *reversed(parts_input))
            assert relative_path.resolve() == absolute_input_path 
            return relative_path
        raise e
    

resources_dir = make_path_relative(meta["resources_dir"])
input1_R1 = resources_dir / "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R1_001.fastq.gz"
input1_R2 = resources_dir / "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_L001_R2_001.fastq.gz"
input2_R1 = resources_dir / "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R1_001.fastq.gz"
input2_R2 = resources_dir / "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_L004_R2_001.fastq.gz"
input3_R1 = resources_dir / "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R1_001.fastq.gz"
input3_R2 = resources_dir / "10x_5k_anticmv/raw/5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_L001_R2_001.fastq.gz"
gex_reference = resources_dir / "reference_gencodev41_chr1/reference_cellranger.tar.gz"
feature_reference = resources_dir / "10x_5k_anticmv/raw/feature_reference.csv"
vdj_reference = resources_dir / "10x_5k_anticmv/raw/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz"

# Beam Input
input1_R1_beam = resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_gex_subset_S3_L001_R1_001.fastq.gz"
input1_R2_beam = resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_gex_subset_S3_L001_R2_001.fastq.gz"
input2_R1_beam = resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_ag_subset_S1_L001_R1_001.fastq.gz"
input2_R2_beam = resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_ag_subset_S1_L001_R2_001.fastq.gz"
input3_R1_beam = resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_vdj_subset_S2_L001_R1_001.fastq.gz"
input3_R2_beam = resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_vdj_subset_S2_L001_R2_001.fastq.gz"
vdj_reference_beam = resources_dir / "10x_5k_beam/raw/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_vdj_reference.tar.gz"
feature_reference_beam =  resources_dir / "10x_5k_beam/raw/beamt_human_A0201_B0702_pbmc_feature_reference.csv"

def test_cellranger_multi(run_component, random_path):
    outputpath = random_path()

    args = [
            "--output", outputpath,
            "--input", input1_R1,
            "--input", input1_R2,
            "--input", input2_R1,
            "--input", input2_R2,
            "--input", input3_R1,
            "--input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
            "--library_type", "Gene Expression;Antibody Capture;VDJ"]
    run_component(args)

    # check for raw data
    assert (outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert (outputpath / "per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()

def test_cellranger_multi_decompressed_reference(run_component, random_path):
    extracted_tar = random_path()
    extracted_tar.mkdir()
    with tarfile.open(gex_reference) as open_tarfile:
        open_tarfile.extractall(extracted_tar)
        run_component([
                "--output", random_path(),
                "--input", input1_R1,
                "--input", input1_R2,
                "--input", input2_R1,
                "--input", input2_R2,
                "--input", input3_R1,
                "--input", input3_R2,
                "--gex_reference", extracted_tar,
                "--vdj_reference", vdj_reference,
                "--feature_reference", feature_reference,
                "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
                "--library_type", "Gene Expression;Antibody Capture;VDJ",
                "--dryrun"])

def test_cellranger_multi_directory_input(run_component, random_path):
    args=[
        "--output", random_path(),
        "--input", meta["resources_dir"] + "/10x_5k_anticmv/raw/",
        "--gex_reference", gex_reference,
        "--vdj_reference", vdj_reference,
        "--feature_reference", feature_reference,
        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
        "--library_type", "Gene Expression;Antibody Capture;VDJ",
        "--gex_secondary_analysis", "true",
        "--gex_generate_bam", "false",
        "--gex_include_introns", "false",
        "--dryrun"]
    run_component(args)

def test_vdj_inner_enrichment_primers(run_component, random_path):
    outputpath = random_path()
    enrichment_primers_file = random_path("txt")
    with enrichment_primers_file.open('w') as primers_file_open:
        primers_file_open.write("AGTCTCTCAGCTGGTACACG\nTCTGATGGCTCAAACACAGC")
    args=[
        "--output", outputpath,
        "--input", meta["resources_dir"] + "/10x_5k_anticmv/raw/",
        "--gex_reference", gex_reference,
        "--vdj_reference", vdj_reference,
        "--feature_reference", feature_reference,
        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
        "--library_type", "Gene Expression;Antibody Capture;VDJ",
        "--gex_secondary_analysis", "true",
        "--gex_generate_bam", "false",
        "--gex_include_introns", "false",
        "--vdj_inner_enrichment_primers", str(make_path_relative(enrichment_primers_file)),
        "--dryrun"]
    run_component(args)
    config_path = outputpath / "config.csv"
    assert config_path.is_file()
    with config_path.open('r') as config_file:
        config_contents = config_file.read()
    expected_csv_content = fr"\[vdj\]\nreference,.*?\ninner-enrichment-primers,{enrichment_primers_file.resolve()}\n"
    assert re.search(expected_csv_content, config_contents)

def test_cellranger_multi_applies_gex_options(run_component, random_path):
    outputpath = random_path()
    args=[
            "--output", outputpath,
            "--input", input1_R1,
            "--input", input1_R2,
            "--input", input2_R1,
            "--input", input2_R2,
            "--input", input3_R1,
            "--input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset;5k_human_antiCMV_T_TBNK_connect_VDJ_subset",
            "--library_type", "Gene Expression;Antibody Capture;VDJ",
            "--gex_secondary_analysis", "true",
            "--gex_generate_bam", "false",
            "--gex_include_introns", "false",
            "--dryrun"]
    run_component(args)
    config_path = outputpath / "config.csv"
    assert config_path.is_file()
    with config_path.open('r') as config_file:
        config_contents = config_file.read()
    expected_csv_content = dedent(
        """\
        chemistry,auto
        no-secondary,False
        create-bam,False
        include-introns,False
        """)
    print (expected_csv_content, flush=True)
    assert expected_csv_content in config_contents

def test_cellranger_multi_no_vdj_reference(run_component, random_path):
    # GH291
    outputpath = random_path()
    args=[
        "--output", outputpath,
        "--input", input1_R1,
        "--input", input1_R2,
        "--input", input2_R1,
        "--input", input2_R2,
        "--input", input3_R1,
        "--input", input3_R2,
        "--gex_reference", gex_reference,
        "--feature_reference", feature_reference,
        "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset;5k_human_antiCMV_T_TBNK_connect_AB_subset",
        "--library_type", "Gene Expression;Antibody Capture",
        "--dryrun"]
    run_component(args)
    assert (outputpath / "config.csv").is_file()

def test_cellranger_multi_crispr_data(run_component, random_path):
    outputpath = random_path()
    args = [
        "--input", meta["resources_dir"] + "/10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_subset_S5_L001_R1_001.fastq.gz",
        "--input", meta["resources_dir"] + "/10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_subset_S5_L001_R2_001.fastq.gz",
        "--input", meta["resources_dir"] + "/10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_subset_S4_L001_R1_001.fastq.gz",
        "--input", meta["resources_dir"] + "/10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_subset_S4_L001_R2_001.fastq.gz",
        "--library_id", "SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_subset;SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_subset",
        "--library_type", "Gene Expression;CRISPR Guide Capture",
        "--min_crispr_umi", "3",
        "--gex_reference", gex_reference,
        "--feature_reference", meta["resources_dir"] + "/10x_5k_lung_crispr/raw/SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference_corrected.csv",
        "--output", outputpath
    ]
    run_component(args)
    # check for raw data
    assert ( outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()
    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()
    # check for filtered gex+ab data
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()
    # check for crispr data
    assert (outputpath / "per_sample_outs/run/count/crispr_analysis/").is_dir()

def test_cellranger_multi_helper_input(run_component, random_path):
    outputpath = random_path()
    args = [
            "--output", outputpath,
            "--gex_input", input1_R1,
            "--gex_input", input1_R2,
            "--abc_input", input2_R1,
            "--abc_input", input2_R2,
            "--vdj_input", input3_R1,
            "--vdj_input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference]
    run_component(args)

    # check for raw data
    assert (outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert (outputpath / "per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()

def test_cellranger_multi_combined_helper_and_global_input(run_component, random_path):
    outputpath = random_path()
    args = [
            "--output", outputpath,
            "--input", input1_R1,
            "--input", input1_R2,
            "--abc_input", input2_R1,
            "--abc_input", input2_R2,
            "--vdj_input", input3_R1,
            "--vdj_input", input3_R2,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset",
            "--library_type", "Gene Expression"]
    run_component(args)

    # check for raw data
    assert (outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert (outputpath / "per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()
    
    
def test_cellranger_multi_create_output_on_fail(run_component, random_path):
    outputpath = random_path()
    # missing vdj_reference
    args = [
            "--output", outputpath,
            "--input", input1_R1,
            "--input", input1_R2,
            "--abc_input", input2_R1,
            "--abc_input", input2_R2,
            "--vdj_input", input3_R1,
            "--vdj_input", input3_R2,
            "--gex_reference", gex_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset",
            "--library_type", "Gene Expression"]
    
    with pytest.raises(subprocess.CalledProcessError) as e:
        run_component(args)
     
    assert (outputpath / "cellranger_multi.log").is_file(), "Should have created log file."


def test_cellranger_multi_beam_data(run_component, random_path):
    outputpath = random_path()
    args = [
        "--input", input1_R1_beam,
        "--input", input1_R2_beam,
        "--input", input2_R1_beam,
        "--input", input2_R2_beam,
        "--input", input3_R1_beam,
        "--input", input3_R2_beam,
        "--library_id", "beamt_human_A0201_B0702_pbmc_gex_subset;beamt_human_A0201_B0702_pbmc_ag_subset;beamt_human_A0201_B0702_pbmc_vdj_subset",
        "--library_type", "Gene Expression;VDJ-T;Antigen Capture",
        "--gex_reference", gex_reference,
        "--vdj_reference", vdj_reference_beam,
        "--feature_reference", feature_reference_beam,
        "--output", outputpath,
        "--control_id", "negative_control_A0201;negative_control_B0702",
        "--mhc_allele", "HLA-A*02:01;HLA-B*07:02"
    ]
    run_component(args)
    # check for raw data
    assert ( outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()
    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()
    # check for antigen data
    assert (outputpath / "per_sample_outs/run/antigen_analysis/").is_dir()
    # check for vdj data
    assert (outputpath / "per_sample_outs/run/vdj_t/").is_dir()


def test_cellranger_multi_fixed_rna(run_component, random_path):
    outputpath = random_path()
    args = [
        "--input", f"{meta['resources_dir']}/10x_5k_fixed/raw/",
        "--library_id", "4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_subset",
        "--library_type", "Gene Expression",
        "--feature_reference", f"{meta['resources_dir']}/10x_5k_fixed/raw/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_feature_reference.csv",
        "--gex_reference", gex_reference,
        "--output", outputpath,
        "--probe_barcode_ids", "BC001;BC002;BC003;BC004",
        "--sample_ids", "Liver_BC1;Ovarian_BC2;Colorectal_BC3;Pancreas_BC4",
        "--gex_generate_bam", "false",
        "--library_lanes", "any",
        "--probe_set", f"{meta['resources_dir']}/10x_5k_fixed/raw/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A_corrected.csv",
        "--sample_force_cells", "5000;-1;-1;-1"
    ]
    run_component(args)
    # check for raw data
    assert (outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()
    # check for metrics summary
    for sample in ["Liver_BC1", "Ovarian_BC2", "Colorectal_BC3", "Pancreas_BC4"]:
        assert (outputpath / f"per_sample_outs/{sample}/metrics_summary.csv").is_file()
        assert (outputpath / f"per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix.h5").is_file()

    assert (outputpath / "multi/multiplexing_analysis").is_dir()


def test_cellranger_multi_with_alternative_names(run_component, random_path):
    import shutil
    import gzip

    input_dir = random_path()
    input_dir.mkdir()

    # Note: if one input file does not use any lanes, none of the input files should use lanes
    # remove lanes
    input1_R1_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_R1_001.fastq.gz"
    input1_R2_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_R2_001.fastq.gz"
    input2_R1_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_R1_001.fastq.gz"
    input2_R2_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_R2_001.fastq.gz"
    input3_R1_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_R1_001.fastq"
    input3_R2_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_R2_001.fastq"

    # copy files
    shutil.copy(input1_R1, input1_R1_link)
    shutil.copy(input1_R2, input1_R2_link)
    shutil.copy(input2_R1, input2_R1_link)
    shutil.copy(input2_R2, input2_R2_link)

    with gzip.open(input3_R1, 'rb') as f_in:
        with open(input3_R1_link, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with gzip.open(input3_R2, 'rb') as f_in:
        with open(input3_R2_link, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    outputpath = random_path()
    args = [
            "--output", outputpath,
            "--input", input1_R1_link,
            "--input", input1_R2_link,
            "--abc_input", input2_R1_link,
            "--abc_input", input2_R2_link,
            "--vdj_input", input3_R1_link,
            "--vdj_input", input3_R2_link,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset",
            "--library_type", "Gene Expression"]
    run_component(args)

    # check for raw data
    assert (outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert (outputpath / "per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()


def test_cellranger_multi_with_alternative_names(run_component, random_path):
    import shutil
    import gzip

    input_dir = random_path()
    input_dir.mkdir()

    # Note: if one input file does not use any lanes, none of the input files should use lanes
    # remove lanes
    input1_R1_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_R1_001.fastq.gz"
    input1_R2_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset_S1_R2_001.fastq.gz"
    input2_R1_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_R1_001.fastq.gz"
    input2_R2_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_AB_subset_S2_R2_001.fastq.gz"
    input3_R1_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_R1_001.fastq"
    input3_R2_link = input_dir / "5k_human_antiCMV_T_TBNK_connect_VDJ_subset_S1_R2_001.fastq"

    # copy files
    shutil.copy(input1_R1, input1_R1_link)
    shutil.copy(input1_R2, input1_R2_link)
    shutil.copy(input2_R1, input2_R1_link)
    shutil.copy(input2_R2, input2_R2_link)

    with gzip.open(input3_R1, 'rb') as f_in:
        with open(input3_R1_link, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    with gzip.open(input3_R2, 'rb') as f_in:
        with open(input3_R2_link, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    outputpath = random_path()
    args = [
            "--output", outputpath,
            "--input", input1_R1_link,
            "--input", input1_R2_link,
            "--abc_input", input2_R1_link,
            "--abc_input", input2_R2_link,
            "--vdj_input", input3_R1_link,
            "--vdj_input", input3_R2_link,
            "--gex_reference", gex_reference,
            "--vdj_reference", vdj_reference,
            "--feature_reference", feature_reference,
            "--library_id", "5k_human_antiCMV_T_TBNK_connect_GEX_1_subset",
            "--library_type", "Gene Expression"]
    run_component(args)

    # check for raw data
    assert (outputpath / "multi/count/raw_feature_bc_matrix.h5").is_file()

    # check for metrics summary
    assert (outputpath / "per_sample_outs/run/metrics_summary.csv").is_file()

    # check for filtered gex+ab data
    assert (outputpath / "per_sample_outs/run/count/sample_filtered_feature_bc_matrix.h5").is_file()

    # check for vdj data
    assert (outputpath / "per_sample_outs/run/vdj_t/filtered_contig_annotations.csv").is_file()

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
