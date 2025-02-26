# openpipelines 2.1.0

## NEW FUNCTIONALITY

* `filter_with_scrublet`: add `expected_doublet_rate`, `stdev_doublet_rate`, `n_neighbors` and `sim_doublet_ratio` arguments (PR #974).

* `feature_annotation/aling_query_reference`: Added a component to align a query and reference dataset (PR #948, #958, #972).

* `workflows/qc/qc` workflow: Added ribosomal gene detection (PR #961).

* `scanvi`: enable CUDA acceleration (PR #969).

* `workflows/annotation/scvi_knn` workflow: Cell-type annotation based on scVI integration followed by KNN label transfer (PR #954).

## MAJOR CHANGES

* Several components: when a component processes a single modality, only that modality is read into memory (PR #944)

* The `transfer/publish` component is deprecated and will be removed in a future major release (PR #941).

# NEW FUNCTIONALITY 

* `workflows/annotation/harmony_knn` workflow: Cell-type annotation based on harmony integration with KNN label transfer (PR #836).

# MINOR CHANGES

* Several workflows: refactor neighbors, leiden and UMAP in a separate subworkflow (PR #942 and PR #949). 

* `grep_annotation_column` and `subset_obsp`: Fix compatibility for SciPy (PR #945).

* `popv`: Pin numpy<2 after new release of scvi-tools (PR #946).

* Various  components (`scgpt` and `annotate`): Add resource labels (PR #947, PR #950).

* `feature_annotation/highly_variable_features_scanpy`: Enable calculation of HVG on a subset of genes (PR #957, PR #959).

* `integrate/scvi`, `integrate/totalvi` and `integrate/scarches`: update base image to nvcr.io/nvidia/pytorch:24.12-py3, pin scvi-tools version to 1.1.5, unpin jax and jaxlib version (PR #970).

* `annotate/celltypist`: Enable passing any layer with log normalized counts, enforce checking whether counts are log normalized (PR #971).

# BUG FIXES

* `utils/subset_vars`: Convert .var column used for subsetting of dtype "boolean" to dtype "bool" when it doesn't contain NaN values (PR #959).

* `resources_test_scripts/annotation_test_data.sh`: Add a layer to the annotation reference dataset with log normalized counts (PR #960).

# openpipelines 2.0.0

## BREAKING CHANGES

* `velocity/scvelo`: update `scvelo` to `0.3.3`, which also removes support for using `loom` input files. The component now uses a `MuData` object as input. Several arguments were added to support selecting different inputs from the MuData file: `counts_layer`, `modality`, `layer_spliced`, `layer_unspliced`, `layer_ambiguous`. An `output_h5mu` argument was has been added (PR #932). 

* `src/annotate/onclass` and `src/annotate/celltypist`: Input parameter for gene name layers of input datasets has been updated to `--input_var_gene_names` and `reference_var_gene_names` (PR #919).

* Several components under `src/scgpt` (`cross_check_genes`, `tokenize_pad`, `binning`) now processes the input (query) datasets differently. Instead of subsetting datasets based on genes in the model vocabulary and/or highly variable genes, these components require an input .var column with a boolean mask specifying this information. The results are written back to the original input data, preserving the dataset structure (PR #832).

* `query/cellxgene_census`: The default output layer has been changed from `.layers["counts"]` to `.X` to be more aligned with the standard OpenPipelines format (PR #933).
  Use argument `--output_layer_counts counts` to revert the behaviour to the previous default.

## NEW FUNCTIONALITY

* `velocyto_to_h5mu`: now writes counts to `.X` (PR #932)

* `qc/calculate_atac_qc_metrics`: new component for calculating ATAC QC metrics (PR #868).

* `workflows/annotation/scgpt_annotation` workflow: Added a scGPT transformer-based cell type annotation workflow (PR #832).

* `workflows/annotation/scgpt_integration_knn` workflow: Cell-type annotation based on scGPT integration with KNN label transfer (PR #875).

* CI: Use `params.resources_test` in test workflows in order to point to an alternative location (e.g. a cache) (PR #889).

## MINOR CHANGES

* Pin `scikit-learn` for `labels_transfer/xgboost` to `<1.6` (PR #931).

* `filter/filter_with_scrublet`: provide cleaner error message when running scrublet on an empty modality (PR #929).

* Several component (cleanup): remove workaround for using being able to use shared utility functions with Nextflow Fusion (PR #920).

* `scgpt/cell_type_annotation` component update: Added support for multi-processing (PR #832).

* Several annotation (`src/annotate/`) components (`onclass`, `celltypist`, `random_forest_annotation`, `scanvi`, `svm_annotation`): Updated input parameteres to ensure uniformity across components, implemented functionality to cross-check the overlap of genes between query and reference (model) datasets and implemented logic to allow for subsetting of genes (PR #919). 

* `workflows/annotation/scgpt_annotation` workflow: Added a scGPT transformer-based cell type annotation workflow (PR #832).

* `scgpt/cross_check_genes` component update: Highly variable genes are now cross-checked based on the boolean mask in `var_input`. The filtering information is stored in the `--output_var_filter` .var field instead of subsetting the dataset (PR #832).

* `scgpt/binning` component update: This component now requires the `--var_input` parameter to provide gene filtering information. Binned data is written to the `--output_obsm_binned_counts` .obsm field in the original input data (PR #832).

* `scgpt/pad_tokenize` component update: Genes are padded and tokenized based on filtering information in `--var_input` and `--input_obsm_binned_counts` (PR #832).

* `resources_test_scripts/scgpt.sh`: Update scGPT test resources to avoid subsetting of datasets (PR #926).

* `workflows/integration/scgpt_leiden` workflow update: Update workflow such that input dataset is not subsetted for HVG but uses boolean masks in .var field instead (PR #875).

## BUG FIXES

* `scvi_leiden` workflow: fix the input layer argument of the workflow not being passed to the scVI component (PR #936 and PR #938). 

* `scgpt/embedding`: remove unused argument `dbsn` (PR #875).

* `scgpt/binning`: update handling of empty rows in sparse matrices (PR #875).

* `dataflow/split_h5mu`: Update memory label from `lowmem` to `highmem` and cpu label from `singlecpu` to `lowcpu` (PR #930).

# openpipelines 2.0.0-rc.2

## BUG FIXES

* `annotate/popv`: fix popv raising `ValueError` when an accelerator (e.g. GPU) is unavailable (PR #915).

## MINOR CHANGES

* `dataflow/split_h5mu`: Optimize resource usage of the component (PR #913).

# openpipelines 2.0.0-rc.1

## BREAKING CHANGES

* Added cell multiplexing support to the `from_cellranger_multi_to_h5mu` component and the `cellranger_multi` workflow. For the `from_cellranger_multi_to_h5mu` component, the `output` argument now requires a value containing a wildcard character `*`, which will be replaced by the sample ID to form the final output file names. Additionally, a `sample_csv` argument is added to the `from_cellragner_multi_to_h5mu` component which describes the sample name per output file. No change is required for the `output_h5mu` argument from the `cellranger_multi` workflow, the workflow will just emit multiple events in case of a multiplexed run, one for each sample. The id of the events (and default output file names) are set by `--sample_ids` (in case of cell multiplexing), or (as before) by the user provided `id` for the input (PR #803 and PR #902).

* `demux/bcl_convert`: update BCL convert from 3.10 to 4.2 (PR #774).

* `demux/cellranger_mkfastq`, `mapping/cellranger_count`, `mapping/cellranger_multi` and `reference/build_cellranger_reference`: update cellranger to `8.0.1` (PR #774 and PR #811).

* Removed `--disable_library_compatibility_check` in favour of `--check_library_compatibility` to the `mapping/cellranger_multi` component and the `ingestion/cellranger_multi` workflow (PR #818).

* `lianapy`: bumped version to `1.3.0` (PR #827 and PR #862). Additionally, `groupby` is now a required argument.

* `concat`: this component was deprecated and has now been removed, use `concatenate_h5mu` instead (PR #796).

* The `workflows` folder in the root of the project no longer contains symbolic links to the build workflows in `target`.
  Using any workflows that was previously linked in this directory will now result in an error which will indicate
  the location of the workflow to be used instead (PR #796).
  
* `XGBoost`: bump version to `2.0.3` (PR #646).

* Several components: update anndata to `0.11.1` and mudata to `0.3.1` (PR #645 and PR #901), and scanpy to `1.10.4` (PR #901). 

* `filter/filter_with_hvg`: this component was deprecated and has now been removed. Use `feature_annotation/highly_variable_features_scanpy` instead (PR #843).

* `dataflow/concat`: this component was deprecated and has now been removed. Use `dataflow/concatenate_h5mu` instead (PR #857).

* `convert/from_h5mu_to_seurat`: bump seurat to latest version (PR #850).

* `workflows/ingestion/bd_rhapsody`: Upgrade BD Rhapsody 1.x to 2.x, thereby changing the interface of the workflow (PR #846).

* `mapping/bd_rhapsody`: Upgrade BD Rhapsody 1.x to 2.x, thereby changing the interface of the workflow (PR #846).

* `reference/make_bdrhap_reference`: Upgrade BD Rhapsody 1.x to 2.x, thereby changing the interface of the workflow (PR #846).

* `reference/build_star_reference`: Rename `mapping/star_build_reference` to `reference/build_star_reference` (PR #846).

* `reference/cellranger_mkgtf`: Rename `reference/mkgtf` to `reference/cellranger_mkgtf` (PR #846).

* `labels_transfer/xgboost`: Align interface with new annotation workflow
  - Store label probabilities instead of uncertainties
  - Take `.h5mu` format as an input instead of `.h5ad`

* `reference/build_cellranger_arc_reference`: a default value of "output" is now specified for the argument `--genome`, inline with `reference/build_cellranger_reference` component. Additionally, providing a value for `--organism` is no longer required and its default value of `Homo Sapiens` has been removed (PR #864).

## MAJOR CHANGES

* Bump popv to `0.4.2` (PR #901)

## NEW FUNCTIONALITY

* Added `demux/cellranger_atac_mkfastq` component: demultiplex raw sequencing data for ATAC experiments (PR #726).

* `process_samples`, `process_batches` and `rna_multisample` workflows: added functionality to scale the log-normalized 
  gene expression data to unit variance and zero mean. The scaled data will be output to a different layer and the
  representation with reduced dimensions will be created and stored in addition to the non-scaled data (PR #733).

* `transform/scaling`: add `--input_layer` and `--output_layer` arguments (PR #733).

* CI: added checking of mudata contents for multiple workflows (PR #783).

* Added multiple arguments to the `cellranger_multi` workflow in order to maintain feature parity with the `mapping/cellranger_multi` component (PR #803).

* `convert/from_cellranger_to_h5mu`: add support for antigen analysis. 

* Added `demux/cellranger_atac_mkfastq` component: demultiplex raw sequencing data for ATAC experiments (PR #726).

* Added `reference/build_cellranger_reference` component: build reference file compatible with ATAC and ATAC+GEX experiments (PR #726).

* `demux/bcl_convert`: add support for no lane splitting (PR #804).

* `reference/cellranger_mkgtf` component: Added cellranger mkgtf as a standalone component (PR #771).

* `scgpt/cross_check_genes` component: Added a gene-model cross check component for scGPT (PR #758).

* `scgpt/embedding`: component: Added scGPT embedding component (PR #761)

* `scgpt/tokenize_pad`: component: Added scGPT padding and tokenization component (PR #754).

* `scgpt/binning` component: Added a scGPT pre-processing binning component (PR #765).

* `workflows/integration/scgpt_leiden` workflow with scGPT integration followed by Leiden clustering (PR #794).

* `scgpt/cell_type_annotation` component: Added scGPT cell type annotation component (PR #798).

* `resources_test_scripts/scGPT.sh`: Added script to include scGPT test resources (PR #800).

* `transform/clr` component: Added the option to set the `axis` along which to apply CLR. Possible to override
  on workflow level as well (PR #767).
  
* `annotate/celltypist` component: Added a CellTypist annotation component (PR #825).

* `dataflow/split_h5mu` component: Added a component to split a single h5mu file into multiple h5mu files based on the values of an .obs column (PR #824).

* `workflows/test_workflows/ingestion` components & `workflows/ingestion`: Added standalone components for integration testing of ingestion workflows (PR #801). 

* `workflows/ingestion/make_reference`: Add additional arguments passed through to the STAR and BD Rhapsody reference components (PR #846).

* `annotate/random_forest_annotation` component: Added a random forest cell type annotation component (PR #848).

* `dataflow/concatenate_h5mu`: data from `.uns`, both originating from the global and per-modality slots, is now retained in the final concatenated output object. Additionally, added the `uns_merge_mode` argument in order to tune the behavior when conflicting keys are detected across samples (PR #859).

* `dimred/densmap` component: Added a densMAP dimensionality reduction component (PR #748).

* `annotate/scanvi` component: Added a component to annotate cells using scANVI (PR #833).

* `transform/bpcells_regress_out` component: Added a component to regress out effects of confounding variables in the count matrix using BPCells (PR #863).

* `transform/regress_out`: Allow providing 'input' and 'output' layers for scanpy regress_out functionality (PR #863).

* `workflows/ingestion/make_reference`: add possibility to build CellRanger ARC references. Added `--motifs_file`, `--non_nuclear_contigs` and `--output_cellranger_arc` arguments (PR #864).

* Test resources (reference_gencodev41_chr1): switch reference genome for CellRanger to ARC variant (PR #864).

* `transform/bpcells_regress_out` component: Added a component to regress out effects of confounding variables in the count matrix using BPCells (PR #863).

* `transform/regress_out`: Allow providing 'input' and 'output' layers for scanpy regress_out functionality (PR #863).

* Added `transform/tfidf` component: normalize ATAC data with TF-IDF (PR #870).

* Added `dimred/lsi` component (PR #552).

* `metadata/duplicate_obs` component: Added a component to make a copy from one .obs field or index to another .obs field within the same MuData object (PR #874, PR #899).

* `annotate/onclass`: component: Added a component to annotate cell types using OnClass (PR #844).

* `annotate/svm` component: Added a component to annotate cell types using support vector machine (SVM) (PR #845).

* `metadata/duplicate_var` component: Added a component to make a copy from one .var field or index to another .var field within the same MuData object (PR #877, PR #899).

* `filter/subset_obsp` component: Added a component to subset an .obsp matrix by column based on the value of an .obs field. The resulting subset is moved to an .obsm field (PR #888).

* `labels_transfer/knn` component: Enable using additional distance functions for KNN classification (PR #830) and allow to perform KNN classification based on a pre-calculated neighborhood graph (PR #890).

## MINOR CHANGES

* Several components: bump python version (PR #901).

* `resources_test_scripts/cellranger_atac_tiny_bcl.sh` script: generate counts from fastq files using CellRanger atac count (PR #726).

* `cellbender_remove_background_v0_2`: update base image to `nvcr.io/nvidia/pytorch:23.12-py3` (PR #646).

* Bump scvelo to `0.3.2` (PR #828).

* Pin numpy<2 for several components (PR #815).

* Added `resources_test_scripts/cellranger_atac_tiny_bcl.sh` script: download tiny bcl file with an ATAC experiment, download a motifs file, demultiplex bcl files to reads in fastq format (PR #726).

* `mapping/cellranger_multi` component now outputs logs on failure of the `cellranger multi` process (PR #766).

* Bump `viash-actions` to `v6` (PR #821).

* `reference/make_reference`: Do not try to extract genome fasta and transcriptome gtf if they are not gzipped (PR #856).

* Changes related to syncing the test resources (PR #867):

  - Add `.info.test_resources` to `_viash.yaml` to specify where test resources need to be synced from.
  - `download/sync_test_resources`: Use `.info.test_resources` in `_viash.yaml` to detect where test resources need to be synced from.
  - Update CI to use `project/sync-and-cache` instead of `project/sync-and-cache-s3`.

## BUG FIXES

* Fix failing tests for `ingestion/cellranger_postprocessing`, `ingestion/conversion` and `multiomics/process_batches` (PR #869).

* `convert/from_10xh5_to_h5mu`: add .uns slot to mdata root when metrics file is provided (PR #887).

* Fix ingestion components not working when optional arguments are unset (PR #894).

* `transform/normalize_total` component: pass the `target_sum` argument to `sc.pp.normalize_total()` (PR #823).

* `from_cellranger_multi_to_h5mu`: fix missing `pytest` dependency (PR #897).

## DOCUMENTATION

* Update authorship of components (PR #835).

# openpipelines 1.0.4

## BUG FIXES

* `scvi_leiden` workflow: fix the input layer argument of the workflow not being passed to the scVI component (PR #939, backported from PR #936 and PR #938). 

# openpipelines 1.0.3

## BUG FIXES

* `qc/calculate_qc_metrics`: increase total counts accuracy with low precision floating dtypes as input layer (PR # , backported from PR #852).

# openpipelines 1.0.2

## BUG FIXES

* `dataflow/concatenate_h5mu`: fix writing out multidimensional annotation dataframes (e.g. `.varm`) that had their 
  data dtype (dtype) changed as a result of adding more observations after concatenation, causing `TypeError`.
  One notable example of this happening is when one of the samples does not have a multimodal annotation dataframe 
  which is present in another sample; causing the values being filled with `NA` (PR #842, backported from PR #837).

# openpipelines 1.0.1

## BUG FIXES

* Bump viash to `0.8.6` (PR #816, backported from #815). This changes the at-runtime generated nextflow process from an in-memory to an on-disk temporary file, which should cause less issues with Nextflow Fusion.

# openpipelines 1.0.0-rc6

## BUG FIXES

* `dataflow/concatenate_h5mu`: fix regression bug where observations are no longer linked to the correct metadata
after concatenation (PR #807)

# openpipelines 1.0.0-rc5

## BUG FIXES

* `cluster/leiden`: prevent leiden component from hanging when a child process is killed (e.g. when there is not enough memory available) (PR #805).

# openpipelines 1.0.0-rc4

## BREAKING CHANGES

* `query/cellxgene_census`: Refactored the interface, documentation and internal workings of this component (PR #621).
  - Renamed arguments to align with standard OpenPipelines notations and cellxgene census API:
    - `--input_database` became `--input_uri`
    - `--cellxgene_release` became `--census_version`
    - `--cell_query` became `--obs_value_filter`
    - `--cells_filter_columns` became `--cell_filter_grouping`
    - `--min_cells_filter_columns` became `--cell_filter_minimum_count`
    - `--modality` became `--output_modality`
    - Removed `--dataset_id` since it was no longer being used.
    - Added `--add_dataset_meta` to add metadata to the output MuData object.
  - Documentation of the component and its arguments was improved.

## BUG FIXES

* `mapping/cellranger_multi`: Fix the regex for the fastq input files to allow dropping the lane from the input file names (e.g. `_L001`) (PR #778).

* `workflows/rna/rna_singlesample`: Fix argument passing `top_n_vars` and `obs_name_mitochondrial_fraction` to the `qc` subworkflow (PR #779).

# openpipelines 1.0.0-rc3

## BREAKING CHANGES

* Docker image names now use `/` instead of `_` between the name of the component and the namespace (PR #712).

## BUG FIXES

* `rna_singlesample`: fixed a bug where selecting the column for the filtering with mitochondrial fractions 
  using `obs_name_mitochondrial_fraction` was done with the wrong column name, causing `ValueError` (PR #743).

* Fix publishing in `process_samples` and `process_batches` (PR #759).

## NEW FUNCTIONALITY

* `dimred/tsne` component: Added a tSNE dimensionality reduction component (PR #742).

# openpipelines 1.0.0-rc2

## BUG FIXES

* Cellranger multi: Fix using a relative input path for `--vdj_inner_enrichment_primers` (PR #717)

* `dataflow/split_modalities`: remove unused `compression` argument. Use `output_compression` instead (PR #714).

* `metadata/grep_annotation_column`: fix calculating fraction when an input observation has no counts, which caused
the result to be out of bounds.

* Fix `--output` argument not working for several workflows (PR #740).

## MINOR CHANGES

* `metadata/grep_annotation_column`: Added more logging output (PR #697).

* `metadata/add_id` and `metadata/grep_annotation_column`: Bump python to 3.11 (PR #697).

* Bump viash to 0.8.5 (PR #697)

* `dataflow/split_modalities`: add more logging output and bump python to 3.12 (PR #714).

* `correction/cellbender`: Update nextflow resource labels from `singlecpu` and `lowmem` to `midcpu` and `midmem` (PR #736)

# openpipelines 1.0.0rc1

## BREAKING CHANGES

* Change separator for arguments with multiple inputs from `:` to `;` (PR #700 and #707). Now, _all_ arguments with `multiple: true` will use `;` as the separator.
  This change was made to be able to deal with file paths that contain `:`, e.g. `s3://my-bucket/my:file.txt`. Furthermore, the `;` separator will become
  the default separator for all arguments with `multiple: true` in Viash >= 0.9.0.

* This project now uses viash version 0.8.4 to build components and workflows. Changes related to this version update should
  be _mostly_ backwards compatible with respect to the results and execution of the pipelines. From a development perspective,
  drastic updates have been made to the developemt workflow.

  Development related changes:
    * Bump viash version to 0.8.4 (PR #598, PR#638 and #706) in the project configuration.
    * All pipelines no longer use the anonymous workflow. Instead, these workflows were given
      a name which was added to the viash config as the entrypoint to the pipeline (PR #598).
    * Removed the `workflows` folder and moved its contents to new locations:
        1. The `resources_test_scripts` folder now resides in the root of the project (PR #605).
        2. All workflows have been moved to the `src/workflows` folder (PR #605).
           This implies that workflows must now be build using `viash (ns) build`, just like with components.
        3. Adjust GitHub Actions to account for new workflow paths (PR #605).
        4. In order to be backwards compatible, the `workflows` folder now contains symbolic
           links to the build workflows in `target`. This is not a problem when using the repository for pipeline
           execution. However, if a developer wishes to contribute to the project, symlink support should be enabled
           in git using `git config core.symlinks=true`. Alternatively, use
           `git clone -c core.symlinks=true git@github.com:openpipelines-bio/openpipeline.git` when cloning the
           repository. This avoids the symlinks being resolved (PR #628). 
        4bis. With PR #668, the workflows have been renamed. This does not hamper the backwards compatibility
              of the symlinks that have been described in 4, because they still use the original location
              which includes the original name.
                * `multiomics/rna_singlesample` has been renamed to `rna/process_single_sample`,
                * `multiomics/rna_multisample` has been renamed to `rna/rna_multisample`,
                * `multiomics/prot_multisample` became `prot/prot_multisample`,
                * `multiomics/prot_singlesample` became `prot/prot_singlesample`,
                * `multiomics/full_pipeline` was moved to `multiomics/process_samples`,
                * `multiomics/multisample` has been renamed to `multiomics/process_batches`,
                * `multiomics/integration/initialize_integration` changed to `multiomics/dimensionality_reduction`,
                * finally, all workflows at `multiomics/integration/*` were moved to `integration/*`

        5. Removed the `workflows/utils` folder. Functionality that was provided by the `DataflowHelper` 
           and `WorkflowHelper` is now being provided by viash when the workflow is being build (PR #605).

  End-user facing changes:
    * The `concat` component had been deprecated and will be removed in a future release.
      It's functionality has been copied to the `concatenate_h5mu` component because the name is in
      conflict with the `concat` operator from nextflow (PR #598).
    * `prot_singlesample`, `rna_singlesample`, `prot_multisample` and `rna_multisample`: QC statistics
       are now only calculated once where needed. This means that the mitochondrial gene detection is
       performed in the `rna_singlesample` pipeline and the other count based statistics are calculated
       during the `prot_multisample` and `rna_multisample` pipelines. In both cases, the `qc` pipeline
       is being used, but only parts of that workflow are activated by parametrization. Previously
       the count based statistics were calculated in both the `singlesample` and `multisample` pipelines,
       with the results from the multisample pipelines overwriting the previous results. What is breaking here
       is that the qc statistics are not being added to the results of the singlesample worklows.
       This is _not_ an issue when using the `full_pipeline` because in this case the singlesample and
       multisample workflows are executed in-tandem. If you wish to execute the singlesample workflows
       in a seperate manner and still include count based statistics, please run the `qc` pipeline
       on the result of the singlesample workflow (PR #604).
    * `filter/filter_with_hvg` has been renamed to `feature_annotation/highly_variable_features_scanpy`, along with the following changes (PR #667).
      - `--do_filter` was removed
      - `--n_top_genes` has been renamed to `--n_top_features`
    * `full_pipeline`, `multisample` and `rna_multisample`: Renamed arguments (PR #667).
      - `--filter_with_hvg_var_output` became `--highly_variable_features_obs_batch_key`
      - `--filter_with_hvg_obs_batch_key` became `--highly_variable_features_var_output`
    * `rna_multisample`: Renamed arguments (PR #667).
      - `--filter_with_hvg_n_top_genes` became `--highly_variable_features_n_top_features`
      - `--filter_with_hvg_flavor` became `--highly_variable_features_flavor`
 
* Renamed `obsm_metrics` to `uns_metrics` for the `cellranger_mapping` workflow because the cellranger metrics are stored in `.uns` and not `.obsm` (PR #610).

## MAJOR CHANGES

* `mapping/cellranger_mkfastq`: update from cellranger `6.0.2` to `7.0.1` (PR #675)

## NEW FUNCTIONALITY

* `multisample` pipeline: This workflow now works when provided multimple unimodal files or multiple multimodal files, in addition to the previously supported single multimodal file (PR #606). The modalities are processed independently from each other:
  - As before, a single multimodal file is split into several unimodal MuData objects, each modality being stored in a file.
  - (New) When multiple unimodal files are provided, they can be used used as is.
  - (New) Mosaic input (i.e. multiple uni- or multimodal files) are split into unimodal files.
    Providing the same modality twice is not supported however, meaning the modalities should be unique.
    For example, using `input: ["data1.h5mu", "data2.h5mu"]` with `data1.h5mu` providing data for `rna` and `atac` 
    and `data2.h5mu` for `rna` and `prot` will not work, because the `rna` modality is present in both input files.
  
* `multisample` workflow: throw an error when argument values for the merge component or the `initialize_integration` workflow differ between the inputs (PR #606).

* Added a `split_modalities` workflow in order to split a multimodal mudata files into several unimodal mudata files. Its behavior is identical to the `split_modalities` component, but it also provides functionality to make sure everything works when nextflow's `-stub` option is enabled (PR #606).

* All workflow now use `dependencies` to handle includes from other workflows (PR #606).

* `qc/calculate_qc_metrics`: allow setting the output column names and disabling the calculation of several metrics (PR #644).

* `rna_multisample`, `prot_multisample` and `qc` workflows: allow setting the output column names and disabling the calculation of several metrics (PR #606).

* `cluster/leiden`: Allow calculating multiple resolutions in parallel (PR #645).

* `qc/calculate_qc_metrics`: allow setting the output column names and disabling the calculation of several metrics (PR #644).

* `rna_multisample` workflow: added `--modality` argument (PR #607).

* `multisample` workflow: in addition to using multimodal files as input, this workflow now also accepts a list of files. The list of files must be the unimodal equivalents of a split multimodal file. The modalities in the list must be unique and after processing the modalities will be merged into multimodal files (PR #606).

* Added `filter/intersect_obs` component which removes observations that are not shared between modalities (PR #589).

* Re-enable `convert/from_h5mu_to_seurat` component (PR #616).

* Added the `gdo_singlesample` pipeline with basic count filtering (PR #672).

* `process_samples` pipeline: the `--rna_layer`, `--prot_layer` and `gdo_layer` argument can not be used to specify an alternative layer to .X where the raw data are stored. To enable this feature, the following changes were required:
  - Added `transform/move_layer` component.
  - `filter/filter_with_scrublet`: added `--layer` argument.
  - `transform/clr`: added `--input_layer` argument.
  - `metadata/grep_annotation_column`: added `--input_layer` argument.
  - `rna/rna_singlesample`, `rna/rna_multisample`, `prot/prot_singlesample` and `prot/prot_multisample`: add `--layer` argument.
  - `process_batches`: Added `rna_layer` and `prot_layer` arguments.

* Enable dataset functionality for nf-tower (PR #701)

* Added `annotate/score_genes` and `annotate/score_genes_cell_cycle` to calculate scanpy gene scores (PR #703).

## MINOR CHANGES

* Refactored `rna_multisample` (PR #607), `cellranger_multi` (PR #609), `cellranger_mapping` (PR #610) and other (PR #606) pipelines to use `fromState` and `toState` functionality.

* `metadata/add_id`: add more runtime logging (PR #663).

* `cluster/leiden`: Bump python to 3.11 and leidenalg to 0.10.0 (PR #645).

* `mapping/htseq_count_to_h5mu` and `multi_star`: update polars and gtfparse (PR #642). 

* Pin `from_h5mu_to_seurat` to use Seurat to version 4 (PR #630).

* `velocity/scvelo`: bump scvelo to 0.3.1 and python to 3.10 (PR #640).

* Updated the Viash YAML schemas to the latest version of Viash (PR #620).

* `build_cellranger_reference` and `build_bdrhap_reference`: Bump go version to `1.21.4` when building seqkit for testing the component (PR #624 and PR #637).

* `correction/cellbender_remove_background`: Remove `muon` as a test dependency (PR #636).

* (Automatic testing) Update viashpy to 0.6.0 (PR #665).

* `integrate/scarches`, `integrate/scvi`, `velocity/scvelo` and `integrate/totalvi`: pin jax, jaxlib to `<0.4.23` (PR #699).

* `integrate/scvi`: Unpin `numba` and pin scvi-tools to `1.0.3` (PR #699).

* `integrate/totalvi`: Enable GPU-accelerated computing, unpin `torchmetrics` and pin jax, jaxlib to `<0.4.23` (PR #699).

## BUG FIXES

* `transform/log1p`: fix `--input_layer` argument not functioning (PR #678). 

* `dataflow/concat` and `dataflow/concatenate_h5mu`: Fix an issue where using `--mode move` on samples with non-overlapping features would cause `var_names` to become unaligned to the data (PR #653).   

* `filter/filter_with_scrublet`: (Testing) Fix duplicate test function names (PR #641).

* `dataflow/concatenate_h5mu` and `dataflow/concat`: Fix `TypeError` when using mode 'move' and a column with conflicting metadata does not exist across all samples (PR #631).

* `dataflow/concatenate_h5mu` and `dataflow/concat`: Fix an issue where joining columns with different datatypes caused `TypeError` (PR #619).

* `qc/calculate_qc_metrics`: Resolved an issue where statistics based on the input columns selected with `--var_qc_metrics` were incorrect when these input columns were encoded in `pd.BooleanDtype()` (PR #685).

* `move_obsm_to_obs`: fix setting output columns when they already exist (PR #690).

* `workflows/dimensionality_reduction` workflow: nearest neighbour calculations no longer recalcalates the PCA when `obm_input` `--obsm_pca` is not set to `X_pca`.

* `feature_annotation/highly_variable_scanpy`: fix .X being used to remove observations with 0 counts when `--layer` has been specified. 

* `filter/filter_with_counts`: fix `--layer` argument not being used.

* `transform/normalize_total`: fix incorrect layer being written to the output when the input layer if not `.X`.

* `src/workflows/qc`: fix input layer not being taken into account when calculating the fraction of mitochondrial genes (always used .X).

* `convert/from_cellranger_multi_to_h5mu`: fix metric values not repesented as percentages being devided by 100. (#704).

# openpipelines 0.12.1

## BUG FIXES

* `rna_singlesample`: Fix filtering parameters values `min_counts`, `max_counts`, `min_genes_per_cell`, `max_genes_per_cell` and `min_cells_per_gene` not being passed to the `filter_with_counts` component (PR #614).

* `prot_singlesample`: Fix filtering parameters values `min_counts`, `max_counts`, `min_proteins_per_cell`, `max_proteins_per_cell` and `min_cells_per_protein` not being passed to the `filter_with_counts` component (PR #614).

# openpipelines 0.12.0

## BREAKING CHANGES

The detection of mitochondrial genes has been revisited in order to remove the interdependency with the count filtering and the QC metric calculation.
Implementing this changes involved breaking some existing functionality:

* `filter/filter_with_counts`: removed `--var_gene_names`, `--mitochondrial_gene_regex`, `--var_name_mitochondrial_genes`, `--min_fraction_mito` and `--max_fraction_mito` (PR #585).

* `workflows/prot_singlesample`: removed `--min_fraction_mito` and `--max_fraction_mito` because regex-based detection detection of mitochondrial genes is not possible (PR #585).

* The fraction of counts that originated from mitochondrial genes used to be written to an .obs column with a name that was derived from `pct_` suffixed by the name of the mitochondrial gene column. The `--obs_name_mitochondrial_fraction` argument is introduced to change the destination column and the default prefix has changed from `pct_` to `fraction_` (PR #585).

## NEW FUNCTIONALITY

* `workflows/qc`: A pipeline to add basic qc statistics to a MuData object (PR #585). 

* `workflows/rna_singlesample`: added `--obs_name_mitochondrial_fraction` and make sure that the values from `--max_fraction_mito`  and `--min_fraction_mito` are bound between 0 and 1 (PR #585).

* Added `filter/delimit_fraction`: Turns an annotation column containing values between 0 and 1 into a boolean column based on thresholds (PR #585).

* Added `metadata/grep_annotation_column`: Perform a regex lookup on a column from the annotation matrices .obs or .var (PR #585).

* `workflows/full_pipelines`: added `--obs_name_mitochondrial_fraction` argument (PR #585).

* `workflows/prot_multisample`: added `--var_qc_metrics` and `--top_n_vars` arguments (PR #585).

* Added genetic demultiplexing methods `cellsnp`, `demuxlet`, `freebayes`, `freemuxlet`, `scsplit`, `sourorcell` and `vireo` (PR #343).

## MINOR CHANGES

* Several components: bump scanpy to 1.9.5 (PR #594).

* Refactored `prot_multisample` and `prot_singlesample` pipelines to use `fromState` and `toState` functionality (PR #585).

# openpipelines 0.11.0

## BREAKING CHANGES

* Nextflow VDSL3: set `simplifyOutput` to `False` by default. This implies that components and workflows will output a hashmap with a sole "output" entry when there is only one output (PR #563).

* `integrate/scvi`: rename `model_output` argument to `output_model` in order to align with the `scvi_leiden` workflow. This also fixes a bug with the workflow where the argument did not function (PR  #562).

## MINOR CHANGES

* `dataflow/concat`: reduce memory consumption when using `--other_axis_mode move` by processing only one annotation matrix (`.var`, `.obs`) at a time (PR #569).

* Update viashpy and pin it to `0.5.0` (PR #572 and PR #577).

* `convert/from_h5ad_to_h5mu`, `convert/from_h5mu_to_h5ad`, `dimred/pca`, `dimred/umap/`, 
`filter/filter_with_counts`, `filter/filter_with_hvg`, `filter/remove_modality`, `filter/subset_h5mu`, 
`integrate/scanorama`, `transform/delete_layer` and `transform/log1p`: update python to `3.9` (PR #572).

* `integrate/scarches`: update base image, `scvi-tools` and `pandas` to `nvcr.io/nvidia/pytorch:23.09-py3`, `~=1.0.3` and `~=2.1.0` respectively (PR #572).

* `integrate/totalvi`: update python to 3.9 and scvi-tools to `~=1.0.3` (PR #572).

* `correction/cellbender_remove_background`: change base image to `nvcr.io/nvidia/cuda:11.8.0-devel-ubuntu22.04` and downwgrade MuData to 0.2.1 because it is the oldest version that uses python 3.7 (PR #575).

* Several integration workflows: prevent leiden from being executed when no resolutions are provided (PR #583).

* `dataflow/concat`: bump pandas to ~=2.1.1 and reduce memory consumption by only reading one modality into memory at a time (PR #568). 

* `annotate/popv`: bump `jax` and `jaxlib` to `0.4.10`, scanpy to `1.9.4`, scvi to `1.0.3` and pin `ml-dtypes` to < 0.3.0 (PR #565).

* `velocity/scvelo`: pin matplotlib to < 3.8.0 (PR #566).

* `mapping/multi_star`: pin multiqc to 1.15.0 (PR #566).

* `mapping/bd_rhapsody`: pin pandas version to <2 (PR #563). 

* `query/cellxgene_census`: replaced label `singlecpu` with label `midcpu`.

* `query/cellxgene_census`: avoid creating MuData object in memory by writing the modality directly to disk (PR #558).

* `integrate/scvi`: use `midcpu` label instead of `singlecpu` (PR #561).

## BUG FIXES

* `transform/clr`: raise an error when CLR fails to return the requested output (PR #579).

* `correction/cellbender_remove_background`: fix missing helper functionality when using Fusion (PR #575).

* `convert/from_bdrhap_to_h5mu`: Avoid `TypeError: Can't implicitly convert non-string objects to strings` by using categorical dtypes when a string column contains NA values (PR #563).

* `qc/calculate_qc_metrics`: fix calculating mitochondrial gene related QC metrics when only or no mitochondrial genes were found (PR #564).

# openpipelines 0.10.1

## MINOR CHANGES

* `integration/scvi_leiden`: Expose hvg selection argument `--var_input` (#543, PR #547).

## BUG FIXES

* `integration/bbknn_leiden`: Set leiden clustering parameter to multiple (#542, PR #545).

* `integration/scvi_leiden`: Fix component name in Viash config (PR #547).

* `integration/harmony_leiden`: Pass `--uns_neighbors` argument `umap` (PR #548).

* Add workaround for bug where resources aren't available when using Nextflow fusion by including `setup_logger`, `subset_vars` and `compress_h5mu` in the script itself (PR #549).

# openpipelines 0.10.0


## BREAKING CHANGES

* `workflows/full_pipeline`: removed `--prot_min_fraction_mito` and `--prot_max_fraction_mito` (PR #451)

* `workflows/rna_multisample` and `workflows/prot_multisample`: Removed concatenation from these pipelines. The input for these pipelines is now a single mudata file that contains data for multiple samples. If you wish to use this pipeline on multiple single-sample mudata files, you can use the `dataflow/concat` components on them first. This also implies that the ability to add ids to multiple single-sample mudata files prior to concatenation is no longer required, hence the removal of `--add_id_to_obs`, `--sample_id`, `--add_id_obs_output`,  and `--add_id_make_observation_keys_unique` (PR #475).

* The `scvi` pipeline was renamed to `scvi_leiden` because `leiden` clustering was added to the pipeline (PR #499).

* Upgrade `correction/cellbender_remove_background` from CellBender v0.2 to CellBender v0.3.0 (PR #523).
  Between these versions, several arguments related to the slots of the output file have been changed.

## MAJOR CHANGES

* Several components: update anndata to 0.9.3 and mudata to 0.2.3 (PR #423).

* Base resources assigned for a process without any labels is now 1 CPU and 2GB (PR #518).

* Updated to Viash 0.7.5 (PR #513).

* Removed deprecated `variant: vdsl3` tags (PR #513).

* Removed unused `version: dev` (PR #513).

* `multiomics/integration/harmony_leiden`: Refactored data flow (PR #513).

* `ingestion/bd_rhapsody`: Refactored data flow (PR #513).

* `query/cellxgene_census`: increased returned metadata content, revised query option, added filtering strategy and refactored functionality (PR #520).

* Refactor loggers using `setup_logger()` helper function (PR #534).

* Refactor unittest tests to pytest tests (PR #534).

## MINOR CHANGES

* Add resource labels to several components (PR #518).

* `full_pipeline`: default value for `--var_qc_metrics` is now the combined values specified for `--mitochondrial_gene_regex` and `--filter_with_hvg_var_output`.

* `dataflow/concat`: reduce memory consumption by only reading one modality at the same time (PR #474).

* Components that use CellRanger, BCL Convert or bcl2fastq: updated from Ubuntu 20.04 to Ubuntu 22.04 (PR #494).

* Components that use CellRanger: updated Picard to 2.27.5 (PR #494).

* `interprete/liana`: Update lianapy to 0.1.9 (PR #497).

* `qc/multiqc`: add unittests (PR #502).

* `reference/build_cellranger_reference`: add unit tests (PR #506).

* `reference/build_bd_rhapsody_reference`: add unittests (PR #504).

## NEW FUNCTIONALITY

* Added `compression/compress_h5mu` component (PR #530).

* Resource management: when a process exits with a status code between 137 and 140, retry the process with increased memory requirements. Memory scales by multiplying the base memory assigned to the process with the attempt number (PR #518 and PR #527).

* `integrate/scvi`: Add `--n_hidden_nodes`, `--n_dimensions_latent_space`, `--n_hidden_layers`, `--dropout_rate`, `--dispersion`, `--gene_likelihood`, `--use_layer_normalization`, `--use_batch_normalization`, `--encode_covariates`, `--deeply_inject_covariates` and `--use_observed_lib_size` parameters.

* `filter/filter_with_counts`: add `--var_name_mitochondrial_genes` argument to store a boolean array corresponding the detected mitochondrial genes.

* `full_pipeline` and `rna_singlesample` pipelines: add `--var_name_mitochondrial_genes`,  `--var_gene_names` and `--mitochondrial_gene_regex` arguments to specify mitochondrial gene detection behaviour.

* `integrate/scvi`: Add `--obs_labels`, `--obs_size_factor`, `--obs_categorical_covariate` and `--obs_continuous_covariate` arguments (PR #496).

* Added `var_qc_metrics_fill_na_value` argument to `calculate_qc_metrics` (PR #477).

* Added `multiomics/multisample` pipeline to run multisample processing followed by the integration setup. It is considered an entrypoint into the full pipeline which skips the single-sample processing. The idea is to allow a a re-run of these steps after a sample has already been processed by the `full_pipeline`. Keep in mind that samples that are provided as input to this pipeline are processed separately and are not concatenated. Hence, the input should be a concatenated sample (PR #475). 

* Added `multiomics/integration/bbknn_leiden` workflow. (PR #456).

* `workflows/prot_multisample` and `workflows/full_pipelines`: add basic QC statistics to prot modality (PR #485).

* `mapping/cellranger_multi`: Add tests for the mapping of Crispr Guide Capture data (PR #494).

* `convert/from_cellranger_multi_to_h5mu`: add `perturbation_efficiencies_by_feature` and `perturbation_efficiencies_by_feature` information to .uns slot of `gdo` modality (PR #494).

* `convert/from_cellranger_multi_to_h5mu`: add `feature_reference` information to the MuData object. Information is split between the modalities. For example `CRISPR Guide Capture` information if added to the `.uns` slot of the `gdo` modality, while `Antibody Capture` information is added to the .uns slot of `prot` (PR #494).

* Added `multiomics/integration/totalvi_leiden` pipeline (PR #500).

* Added totalVI component (PR #386).

* `workflows/full_pipeline`: Add `pca_overwrite` argument (PR #511).

* Add `main_build_viash_hub` action to build, tag, and push components and docker images for viash-hub.com (PR #480).

* `integration/bbknn_leiden`: Update state management to `fromState` / `toState` (PR #538).

* `mapping/cellranger_multi`: Add optional helper input: allow for passing modality specific inputs, from which library type and library id are inferred (PR #693).

## DOCUMENTATION

* `images`: Added images for various concepts, such as a sample, a cell, RNA, ADT, ATAC, VDJ (PR #515).

* `multiomics/rna_singlesample`: Add image for workflow (PR #515).

* `multiomics/rna_multisample`: Add image for workflow (PR #515).

* `multiomics/prot_singlesample`: Add image for workflow (PR #515).

* `multiomics/prot_multisample`: Add image for workflow (PR #515).

## BUG FIXES

* Fix an issue with `workflows/multiomics/scanorama_leiden` where the `--output` argument doesn't work as expected (PR #509).

* Fix an issue with `workflows/full_pipeline` not correctly caching previous runs (PR #460).

* Fix incorrect namespaces of the integration pipelines (PR #464).

* Fix an issue in several workflows where the `--output` argument would not work (PR #476).

* `integration/harmony_leiden` and `integration/scanorama_leiden`: Fix an issue where the prefix of the columns that store the leiden clusters was hardcoded to `leiden`, instead of adapting to the value for `--obs_cluster` (PR #482). 

* `velocity/velocyto`: Resolve symbolic link before checking whether the transcriptome is a gzip (PR #484).

* `workflows/integration/scanorama_leiden`: fix an issue where `--obsm_input`, --obs_batch`, `--batch_size`, `--sigma`, `--approx`, `--alpha` and `-knn` were not working beacuse they were not passed through to the scanorama component (PR #487).

* `workflows/integration/scanorama_leiden`: fix leiden being calculated on the wrong embedding because the `--obsm_input` argument was not correctly set to the output embedding of scanorama (PR #487).

* `mapping/cellranger_multi`: Fix and issue where modalities did not have the proper name (PR #494).

* `metadata/add_uns_to_obs`: Fix `KeyError: 'ouput_compression'` error (PR #501).

* `neighbors/bbknn`: Fix `--input` not being a required argument (PR #518).

* Create `correction/cellbender_remove_background_v0.2` for legacy CellBender v0.2 format (PR #523).

* `integrate/scvi`: Ensure output has the same dimensionality as the input (PR #524).

* `mapping/bd_rhapsody`: Fix `--dryrun` argument not working (PR #534).

* `qc/multiqc`: Fix component not working for multiple inputs (PR #537). Also converted Bash script to Python scripts.

* `neighbors/bbknn`: Fix `--uns_output`, `--obsp_distances` and `--obsp_connectivities` not being processed correctly (PR #538).

# openpipelines 0.9.0

## BREAKING CHANGES

Running the integration in the `full_pipeline` deemed to be impractical because a plethora of integration methods exist, which in turn interact with other functionality (like clustering). This generates a large number of possible usecases which one pipeline cannot cover in an easy manner. Instead, each integration methods will be split into its separate pipeline, and the `full_pipeline` will prepare for integration by performing steps that are required by many integration methods. Therefore, the following changes were performed:

  * `workflows/full_pipeline`: `harmony` integration and `leiden` clustering are removed from the pipeline.

  * Added `initialize_integration` to run calculations that output information commonly required by the integration methods. This pipeline runs PCA, nearest neighbours and UMAP. This pipeline is run as a subpipeline at the end of `full_pipeline`.

  * Added `leiden_harmony` integration pipeline: run harmony integration followed by neighbour calculations and leiden clustering. Also runs umap on the result.

  * Removed the `integration` pipeline.

The old behavior of the `full_pipeline` can be obtained by running `full_pipeline` followed by the `leiden_harmony` pipeline.

* The `crispr` and `hashing` modalities have been renamed to `gdo` and `hto` respectively (PR #392).

* Updated Viash to 0.7.4 (PR #390).

* `cluster/leiden`: Output is now stored into `.obsm` instead of `.obs` (PR #431).

## NEW FUNCTIONALITY

* `cluster/leiden` and `integration/harmony_leiden`: allow running leiden multiple times with multiple resolutions (PR #431).

* `workflows/full_pipeline`: PCA, nearest neighbours and UMAP are now calculated for the `prot` modality (PR #396).

* `transform/clr`: added `output_layer` argument (PR #396).

* `workflows/integration/scvi`: Run scvi integration followed by neighbour calculations and run umap on the result (PR #396).

* `mapping/cellranger_multi` and `workflows/ingestion/cellranger_multi`: Added `--vdj_inner_enrichment_primers` argument (PR #417).

* `metadata/move_obsm_to_obs`: Move a matrix from an `.obsm` slot into `.obs` (PR #431).

* `integrate/scvi` validity checks for non-normalized input, obs and vars in order to proceed to training (PR #429).

* `schemas`: Added schema files for authors (PR #436).

* `schemas`: Added schema file for Viash configs (PR #436).

* `schemas`: Refactor author import paths (PR #436).

* `schemas`: Added schema file for file format specification files (PR #437).

* `query/cellxgene_census`: Query Cellxgene census component and save the results to a MuData file. (PR #433).

## MAJOR CHANGES

* `report/mermaid`: Now used `mermaid-cli` to generate images instead of creating a request to `mermaid.ink`. New `--output_format`, `--width`, `--height` and  `--background_color` arguments were added (PR #419).

* All components that used `python` as base container: use `slim` version to reduce container image size (PR #427).

## MINOR CHANGES

* `integrate/scvi`: update scvi to 1.0.0 (PR #448)

* `mapping/multi_star`: Added `--min_success_rate` which causes component to fail when the success rate of processed samples were successful (PR #408).

* `correction/cellbender_remove_background` and `transform/clr`: update muon to 0.1.5 (PR #428)

* `ingestion/cellranger_postprocessing`: split integration tests into several workflows (PR #425).

* `schemas`: Add schema file for author yamls (PR #436).

* `mapping/multi_star`, `mapping/star_build_reference` and `mapping/star_align`: update STAR from 2.7.10a to 2.7.10b (PR #441).

## BUG FIXES

* `annotate/popv`: Fix concat issue when the input data has multiple layers (#395, PR #397).

* `annotate/popv`: Fix indexing issue when MuData object contain non overlapping modalities (PR #405).

* `mapping/multi_star`: Fix issue where temp dir could not be created when group_id contains slashes (PR #406).

* `mapping/multi_star_to_h5mu`: Use glob to look for count files recursively (PR #408).

* `annotate/popv`: Pin `PopV`, `jax` and `jaxlib` versions (PR #415).

* `integrate/scvi`: the max_epochs is no longer required since it has a default value (PR #396).

* `workflows/full_pipeline`: fix `make_observation_keys_unique` parameter not being correctly passed to the `add_id` component, causing `ValueError: Observations are not unique across samples` during execution of the `concat` component (PR #422).

* `annotate/popv`: now sets `aprox` to `False` to avoid using `annoy` in scanorama because it fails on processors that are missing the AVX-512 instruction sets, causing `Illegal instruction (core dumped)`.

* `workflows/full_pipeline`: Avoid adding sample names to observation ids twice (PR #457). 

# openpipelines 0.8.0

## BREAKING CHANGES

* `workflows/full_pipeline`: Renamed inconsistencies in argument naming (#372):
  - `rna_min_vars_per_cell` was renamed to `rna_min_genes_per_cell`
  - `rna_max_vars_per_cell` was renamed to `rna_max_genes_per_cell`
  - `prot_min_vars_per_cell` was renamed to `prot_min_proteins_per_cell`
  - `prot_max_vars_per_cell` was renamed to `prot_max_proteins_per_cell`

* `velocity/scvelo`: bump anndata from <0.8 to 0.9.

## NEW FUNCTIONALITY

* Added an extra label `veryhighmem` mostly for `cellranger_multi` with a large number of samples.

* Added `multiomics/prot_multisample` pipeline.

* Added `clr` functionality to `prot_multisample` pipeline.

* Added `interpret/lianapy`: Enables the use of any combination of ligand-receptor methods and resources, and their consensus.

* `filter/filter_with_scrublet`: Add `--allow_automatic_threshold_detection_fail`: when scrublet fails to detect doublets, the component will now put `NA` in the output columns.

* `workflows/full_pipeline`: Allow not setting the sample ID to the .obs column of the MuData object.

* `workflows/rna_multisample`: Add the ID of the sample to the .obs column of the MuData object.

* `correction/cellbender_remove_background`: add `obsm_latent_gene_encoding` parameter to store the latent gene representation.

## BUG FIXES

* `transform/clr`: fix anndata object instead of matrix being stored as a layer in output `MuData`, resulting in `NoneTypeError` object after reading the `.layers` back in.

* `dataflow/concat` and `dataflow/merge`: fixed a bug where boolean values were cast to their string representation.

* `workflows/full_pipeline`: fix running pipeline with `-stub`.

* Fixed an issue where passing a remote file URI (for example `http://` or `s3://`) as `param_list` caused `No such file` errors.

* `workflows/full_pipeline`: Fix incorrectly named filtering arguments (#372).

* `integrate/scvi`: Fix bug when subsetting using the `var_input` argument (PR #385).
* 
* `correction/cellbender_remove_background`: add `obsm_latent_gene_encoding` parameter to store the latent gene representation.

## MINOR CHANGES

* `integrate/scarches`, `integrate/scvi` and `correction/cellbender_remove_background`: Update base container to `nvcr.io/nvidia/pytorch:22.12-py3`

* `integrate/scvi`: add `gpu` label for nextflow platform.

* `integrate/scvi`: use cuda enabled `jax` install.

* `convert/from_cellranger_multi_to_h5mu`, `dataflow/concat` and `dataflow/merge`: update pandas to 2.0.0

* `dataflow/concat` and `dataflow/merge`: Boolean and integer columns are now represented by the `BooleanArray` and `IntegerArray` dtypes in order to allow storing `NA` values.

* `interpret/lianapy`: use the latest development release (commit 11156ddd0139a49dfebdd08ac230f0ebf008b7f8) of lianapy in order to fix compatibility with numpy 1.24.x.

* `filter/filter_with_hvg`: Add error when specified input layer cannot be found in input data.

* `workflows/multiomics/full_pipeline`: publish the output from sample merging to allow running different integrations.

* CI: Remove various unused software libraries from runner image in order to avoid `no space left on device` (PR #425, PR #447).

# openpipelines 0.7.1

## NEW FUNCTIONALITY

* `integrate/scvi`: use `nvcr.io/nvidia/pytorch:22.09-py3` as base container to enable GPU acceleration.

* `integrate/scvi`: add `--model_output` to save model.

* `workflows/ingestion/cellranger_mapping`: Added `output_type` to output the filtered Cell Ranger data as h5mu, not the converted raw 10xh5 output.

* Several components:  added `--output_compression` component to set the compression of output .h5mu files.

* `workflows/full_pipeline` and `workflows/integration`: Added `leiden_resolution` argument to control the coarseness of the clustering.

* Added `--rna_theta` and `--rna_harmony_theta` to full and integration pipeline respectively in order to tune the diversity clustering penalty parameter for harmony integration.

* `dimred/pca`: fix `variance` slot containing a second copy of the variance ratio matrix and not the variances.

## BUG FIXES

* `mapping/cellranger_multi`: Fix an issue where using a directory as value for `--input` would cause `AttributeError`.

* `workflows/integration`: `init_pos` is no longer set to the integration layer (e.g. `X_pca_integrated`).

## MINOR CHANGES

* `integration` and `full` workflows: do not run harmony integration when `obs_covariates` is not provided.

* Add `highmem` label to `dimred/pca` component.

* Remove disabled `convert/from_csv_to_h5mu` component.

* Update to Viash 0.7.1.

* Several components: update to scanpy 1.9.2

* `process_10xh5/filter_10xh5`: speed up build by using `eddelbuettel/r2u:22.04` base container.

## MAJOR CHANGES

* `dataflow/concat`: Renamed `--compression` to `--output_compression`.

# openpipelines 0.7.0

## MAJOR CHANGES

* Removed `bin` folder. As of viash 0.6.4, a `_viash.yaml` file can be included in the root of a repository to set common viash options for the project.
These options were previously covered in the `bin/init` script, but this new feature of viash makes its use unnecessary. The `viash` and `nextlow` should now be installed in a directory that is included in your `$PATH`.

## MINOR CHANGES

* `filter/do_filter`: raise an error instead of printing a warning when providing a column for `var_filer` or `obs_filter` that doesn't exist.

## BUG FIXES

* `workflows/full_pipeline`: Fix setting .var output column for filter_with_hvg.

* Fix running `mapping/cellranger_multi` without passing all references.

* `filter/filter_with_scrublet`: now sets `use_approx_neighbors` to `False` to avoid using `annoy` because it fails on processors that are missing the AVX-512 instruction sets.

* `workflows`: Updated `WorkflowHelper` to newer version that allows applying defaults when calling a subworkflow from another workflow.

* Several components: pin matplotlib to <3.7 to fix scanpy compatibility (see https://github.com/scverse/scanpy/issues/2411).  

* `workflows`: fix a bug when running a subworkflow from a workflow would cause the parent config to be read instead of the subworklow config.

* `correction/cellbender_remove_background`: Fix description of input for cellbender_remove_background.

* `filter/do_filter`: resolved an issue where the .obs column instead of the .var column was being logged when filtering using the .var column.

* `workflows/rna_singlesample` and `workflows/prot_singlesample`: Correctly set var and obs columns while filtering with counts.

* `filter/do_filter`: removed the default input value for `var_filter` argument.

* `workflows/full_pipeline` and `workflows/integration`: fix PCA not using highly variable genes filter.

# openpipelines 0.6.2

## NEW FUNCTIONALITY

* `workflows/full_pipeline`: added `filter_with_hvg_obs_batch_key` argument for batched detection of highly variable genes.

* `workflows/rna_multisample`: added `filter_with_hvg_obs_batch_key`, `filter_with_hvg_flavor` and `filter_with_hvg_n_top_genes` arguments.

* `qc/calculate_qc_metrics`: Add basic statistics: `pct_dropout`, `num_zero_obs`, `obs_mean` and `total_counts` are added to .var. `num_nonzero_vars`, `pct_{var_qc_metrics}`, `total_counts_{var_qc_metrics}`, `pct_of_counts_in_top_{top_n_vars}_vars` and `total_counts` are included in .obs

* `workflows/multiomics/rna_multisample` and `workflows/multiomics/full_pipeline`: add `qc/calculate_qc_metrics` component to workflow.

* `workflows/multiomics/prot_singlesample`: Processing unimodal single-sample CITE-seq data.

* `workflows/multiomics/rna_singlesample` and `workflows/multiomics/full_pipeline`: Add filtering arguments to pipeline.

## MINOR CHANGES

* `convert/from_bdrhap_to_h5mu`: bump R version to 4.2.

* `process_10xh5/filter_10xh5`: bump R version to 4.2.

* `dataflow/concat`: include path of file in error message when reading a mudata file fails.

* `mapping/cellranger_multi`: write cellranger console output to a `cellranger_multi.log` file.

## BUG FIXES

* `mapping/htseq_count_to_h5mu`: Fix a bug where reading in the gtf file caused `AttributeError`. 

* `dataflow/concat`: the `--input_id` is no longer required when `--mode` is not `move`.

* `filter/filter_with_hvg`: does no longer try to use `--varm_name` to set non-existant metadata when running with `--flavor seurat_v3`, which was causing `KeyError`.

* `filter/filter_with_hvg`: Enforce that `n_top_genes` is set when `flavor` is set to 'seurat_v3'.

* `filter/filter_with_hvg`: Improve error message when trying to use 'cell_ranger' as `flavor` and passing unfiltered data.

* `mapping/cellranger_multi` now applies `gex_chemistry`, `gex_secondary_analysis`, `gex_generate_bam`, `gex_include_introns` and `gex_expect_cells`.

# openpipeline 0.6.1

## NEW FUNCTIONALITY

* `mapping/multi_star`: A parallellized version of running STAR (and HTSeq).

* `mapping/multi_star_to_h5mu`: Convert the output of `multi_star` to a h5mu file.

## BUG FIXES

* `filter/filter_with_counts`: Fix an issue where mitochrondrial genes were being detected in .var_names, which contain ENSAMBL IDs instead of gene symbols in the pipelines. Solution was to create a `--var_gene_names` argument which allows selecting a .var column to check using a regex (`--mitochondrial_gene_regex`).

* `dataflow/concat`, `report/mermaid`, `transform/clr`: Don't forget to exit with code returned by pytest.
# openpipeline 0.6.0

## NEW FUNCTIONALITY

* `workflows/full_pipeline`: add `filter_with_hvg_var_output` argument.

* `dimred/pca`: Add `--overwrite` and `--var_input` arguments.

* `tranform/clr`: Perform CLR normalization on CITE-seq data.

* `workflows/ingestion/cellranger_multi`: Run Cell Ranger multi and convert the output to .h5mu.

* `filter/remove_modality`: Remove a single modality from a MuData file.

* `mapping/star_align`: Align `.fastq` files using STAR.

* `mapping/star_align_v273a`: Align `.fastq` files using STAR v2.7.3a.

* `mapping/star_build_reference`: Create a STAR reference index.

* `mapping/cellranger_multi`: Align fastq files using Cell Ranger multi.

* `mapping/samtools_sort`: Sort and (optionally) index alignments.

* `mapping/htseq_count`: Quantify gene expression for subsequent testing for differential expression.

* `mapping/htseq_count_to_h5mu`: Convert one or more HTSeq outputs to a MuData file.

* Added from `convert/from_cellranger_multi_to_h5mu` component.

## MAJOR CHANGES

* `convert/from_velocyto_to_h5mu`: Moved to `velocity/velocyto_to_h5mu`.
  It also now accepts an optional `--input_h5mu` argument, to allow directly reading
  the RNA velocity data into a `.h5mu` file containing the other modalities.

* `resources_test/cellranger_tiny_fastq`: Include RNA velocity computations as part of
  the script.

* `mapping/cellranger_mkfastq`: remove --memory and --cpu arguments as (resource management is automatically provided by viash).

## MINOR CHANGES

* Several components: use `gzip` compression for writing .h5mu files.

* Default value for `obs_covariates` argument of full pipeline is now `sample_id`.

* Set the `tag` directive of all Nextflow components to '$id'.

## BUG FIXES

* Keep data for modalities that are not specifically enabled when running full pipeline.

* Fix many components thanks to Viash 0.6.4, which causes errors to be 
  thrown when input and output files are defined but not found.


# openpipeline 0.5.1

## BREAKING CHANGES

* `reference/make_reference`: Input files changed from `type: string` to `type: file` to allow Nextflow to cache the input files fetched from URL.

* several components (except `from_h5ad_to_h5mu`): the `--modality` arguments no longer accept multiple values.

* Remove outdated `resources_test_scripts`.

* `convert/from_h5mu_to_seurat`: Disabled because MuDataSeurat is currently broken, see [https://github.com/PMBio/MuDataSeurat/issues/9](PMBio/MuDataSeurat#9).

* `integrate/harmony`: Disabled because it is currently not functioning and the alternative, harmonypy, is used in the workflows.

* `dataflow/concat`: Renamed --sample_names to --input_id and moved the ability to add sample id and to join the sample ids with the observation names to `metadata/add_id`

* Moved `dataflow/concat`, `dataflow/merge` and `dataflow/split_modalities` to a new namespace: `dataflow`.

* Moved `workflows/conversion/conversion` to `workflows/ingestion/conversion`

## NEW FUNCTIONALITY

* `metadata/add_id`: Add an id to a column in .obs. Also allows joining the id to the .obs_names.

* `workflows/ingestion/make_reference`: A generic component to build a transcriptomics reference into one of many formats.

* `integrate/scvi`: Performs scvi integration.

* `integrate/add_metadata`: Add a csv containing metadata to the .obs or .var field of a mudata file.

* `DataflowHelper.nf`: Added `passthroughMap`. Usage:

  ```groovy
  include { passthroughMap as pmap } from "./DataflowHelper.nf"
  
  workflow {
    Channel.fromList([["id", [input: "foo"], "passthrough"]])
      | pmap{ id, data ->
        [id, data + [arg: 10]]
      }
  }
  ```
  Note that in the example above, using a regular `map` would result in an exception being thrown,
  that is, "Invalid method invocation `call` with arguments".

  A synonymous of doing this with a regular `map()` would be:
  ```groovy
  workflow {
    Channel.fromList([["id", [input: "foo"], "passthrough"]])
      | map{ tup ->
        def (id, data) = tup
        [id, data + [arg: 10]] + tup.drop(2)
      }
  }
  ```

* `correction/cellbender_remove_background`: Eliminating technical artifacts from high-throughput single-cell RNA sequencing data.

* `workflows/ingestion/cellranger_postprocessing`: Add post-processing of h5mu files created from Cell Ranger data.

* `annotate/popv`: Performs popular major vote cell typing on single cell sequence data.

## MAJOR CHANGES

* `workflows/utils/DataflowHelper.nf`: Added helper functions `setWorkflowArguments()` and `getWorkflowArguments()` to split the data field of a channel event into a hashmap. Example usage:
  ```groovy
  | setWorkflowArguments(
    pca: [ "input": "input", "obsm_output": "obsm_pca" ]
    integration: [ "obs_covariates": "obs_covariates", "obsm_input": "obsm_pca" ]
  )
  | getWorkflowArguments("pca")
  | pca
  | getWorkflowArguments("integration")
  | integration
  ```

* `mapping/cellranger_count`: Allow passing both directories as well as individual fastq.gz files as inputs.

* `convert/from_10xh5_to_h5mu`: Allow reading in QC metrics, use gene ids as `.obs_names` instead of gene symbols.

* `workflows/conversion`: Update pipeline to use the latest practices and to get it to a working state.

## MINOR CHANGES

* `dimred/umap`: Streamline UMAP parameters by adding `--obsm_output` parameter to allow choosing the output `.obsm` slot.

* `workflows/multiomics/integration`: Added arguments for tuning the various output slots of the integration pipeline, namely `--obsm_pca`, `--obsm_integrated`, `--uns_neighbors`, `--obsp_neighbor_distances`, `--obsp_neighbor_connectivities`, `--obs_cluster`, `--obsm_umap`.

* Switch to Viash 0.6.1.

* `filter/subset_h5mu`: Add `--modality` argument, export to VDSL3, add unit test.

* `dataflow/split_modalities`: Also output modality types in a separate csv.

## BUG FIXES

* `convert/from_bd_to_10x_molecular_barcode_tags`: Replaced UTF8 characters with ASCII. OpenJDK 17 or lower might throw the following exception when trying to read a UTF8 file: `java.nio.charset.MalformedInputException: Input length = 1`.

* `dataflow/concat`: Overriding sample name in .obs no longer raises `AttributeError`.

* `dataflow/concat`: Fix false positives when checking for conflicts in .obs and .var when using `--mode move`.

# openpipeline 0.5.0

Major redesign of the integration and multiomic workflows. Current list of workflows:

* `ingestion/bd_rhapsody`: A generic pipeline for running BD Rhapsody WTA or Targeted mapping, with support for AbSeq, VDJ and/or SMK.

* `ingestion/cellranger_mapping`: A pipeline for running Cell Ranger mapping.

* `ingestion/demux`: A generic pipeline for running bcl2fastq, bcl-convert or Cell Ranger mkfastq.

* `multiomics/rna_singlesample`: Processing unimodal single-sample RNA transcriptomics data.

* `multiomics/rna_multisample`: Processing unimodal multi-sample RNA transcriptomics data.

* `multiomics/integration`: A pipeline for demultiplexing multimodal multi-sample RNA transcriptomics data.

* `multiomics/full_pipeline`: A pipeline to analyse multiple multiomics samples.

## BREAKING CHANGES

* Many components: Renamed `.var["gene_ids"]` and `.var["feature_types"]` to `.var["gene_id"]` and `.var["feature_type"]`.

## DEPRECATED

* `convert/from_10xh5_to_h5ad` and `convert/from_bdrhap_to_h5ad`: Removed h5ad based components.

* `mapping/bd_rhapsody_wta` and `workflows/ingestion/bd_rhapsody_wta`: Deprecated in favour for more generic `mapping/bd_rhapsody` and `workflows/ingestion/bd_rhapsody` pipelines.

* `convert/from_csv_to_h5mu`: Disable until it is needed again.

* `dataflow/concat`: Deprecated `"concat"` option for `--other_axis_mode`.

## NEW COMPONENTS

* `graph/bbknn`: Batch balanced KNN.

* `transform/scaling`: Scale data to unit variance and zero mean.

* `mapping/bd_rhapsody`: Added generic component for running the BD Rhapsody WTA or Targeted analysis, with support for AbSeq, VDJ and/or SMK.

* `integrate/harmony` and `integrate/harmonypy`: Run a Harmony integration analysis (R-based and Python-based, respectively).

* `integrate/scanorama`: Use Scanorama to integrate different experiments.

* `reference/make_reference`: Download a transcriptomics reference and preprocess it (adding ERCC spikeins and filtering with a regex).

* `reference/build_bdrhap_reference`: Compile a reference into a STAR index in the format expected by BD Rhapsody.

## NEW WORKFLOWS

* `workflows/ingestion/bd_rhapsody`: Added generic workflow for running the BD Rhapsody WTA or Targeted analysis, with support for AbSeq, VDJ and/or SMK.

* `workflows/multiomics/full_pipeline`: Implement pipeline for processing multiple multiomics samples.

## NEW FUNCTIONALITY

* `convert/from_bdrhap_to_h5mu`: Added support for being able to deal with WTA, Targeted, SMK, AbSeq and VDJ data.

* `dataflow/concat`: Added `"move"` option to `--other_axis_mode`, which allows merging `.obs` and `.var` by only keeping elements of the matrices which are the same in each of the samples, moving the conflicting values to `.varm` or `.obsm`.

## MAJOR CHANGES

* Multiple components: Update to anndata 0.8 with mudata 0.2.0. This means that the format of the `.h5mu` files have changed.

* `multiomics/rna_singlesample`: Move transformation counts into layers instead of overwriting `.X`.

* Updated to Viash 0.6.0.

## MINOR CHANGES

* `velocity/velocyto`: Allow configuring memory and parallellisation.

* `cluster/leiden`: Add `--obsp_connectivities` parameter to allow choosing the output slot.

* `workflows/multiomics/rna_singlesample`, `workflows/multiomics/rna_multisample` and `workflows/multiomics/integration`: Allow choosing the output paths.

* `neighbors/bbknn` and `neighbors/find_neighbors`: Add parameters for choosing the input/output slots.

* `dimred/pca` and `dimred/umap`: Add parameters for choosing the input/output slots.

* `dataflow/concat`: Optimize concat performance by adding multiprocessing and refactoring functions.

* `workflows/multimodal_integration`: Add `obs_covariates` argument to pipeline.

## BUG FIXES

* Several components: Revert using slim versions of containers because they do not provide the tools to run nextflow with trace capabilities.

* `dataflow/concat`: Fix an issue where joining boolean values caused `TypeError`.

* `workflows/multiomics/rna_multisample`, `workflows/multiomics/rna_singlesample` and `workflows/multiomics/integration`: Use nextflow trace reporting when running integration tests.


# openpipeline 0.4.1

## BUG FIXES

* `workflows/ingestion/bd_rhapsody_wta`: use ':' as a seperator for multiple input files and fix integration tests.

## MINOR CHANGES

* Several components: pin mudata and scanpy dependencies so that anndata version <0.8.0 is used.

# openpipeline 0.4.0

## NEW FUNCTIONALITY

* `convert/from_bdrhap_to_h5mu`: Merge one or more BD rhapsody outputs into an h5mu file.

* `dataflow/split_modalities`: Split the modalities from a single .h5mu multimodal sample into seperate .h5mu files. 

* `dataflow/concat`: Combine data from multiple samples together.

## MINOR CHANGES

* `mapping/bd_rhapsody_wta`: Update to BD Rhapsody 1.10.1.

* `mapping/bd_rhapsody_wta`: Add parameters for overriding the minimum RAM & cores. Add `--dryrun` parameter.

* Switch to Viash 0.5.14.

* `convert/from_bdrhap_to_h5mu`: Update to BD Rhapsody 1.10.1.

* `resources_test/bdrhap_5kjrt`: Add subsampled BD rhapsody datasets to test pipeline with.

* `resources_test/bdrhap_ref_gencodev40_chr1`: Add subsampled reference to test BD rhapsody pipeline with.

* `dataflow/merge`: Merge several unimodal .h5mu files into one multimodal .h5mu file.

* Updated several python docker images to slim version.

* `mapping/cellranger_count_split`: update container from ubuntu focal to ubuntu jammy

* `download/sync_test_resources`: update AWS cli tools from 2.7.11 to 2.7.12 by updating docker image

* `download/download_file`: now uses bash container instead of python.

* `mapping/bd_rhapsody_wta`: Use squashed docker image in which log4j issues are resolved.

## BUG FIXES

* `workflows/utils/WorkflowHelper.nf`: Renamed `utils.nf` to `WorkflowHelper.nf`.

* `workflows/utils/WorkflowHelper.nf`: Fix error message when required parameter is not specified.

* `workflows/utils/WorkflowHelper.nf`: Added helper functions:
  - `readConfig`: Read a Viash config from a yaml file.
  - `viashChannel`: Create a channel from the Viash config and the params object.
  - `helpMessage`: Print a help message and exit.

* `mapping/bd_rhapsody_wta`: Update picard to 2.27.3.

## DEPRECATED

* `convert/from_bdrhap_to_h5ad`: Deprecated in favour for `convert/from_bdrhap_to_h5mu`.

* `convert/from_10xh5_to_h5ad`: Deprecated in favour for `convert/from_10xh5_to_h5mu`.

# openpipeline 0.3.1

## NEW FUNCTIONALITY

* `bin/port_from_czbiohub_utilities.sh`: Added helper script to import components and pipelines from `czbiohub/utilities`

Imported components from `czbiohub/utilities`:

* `demux/cellranger_mkfastq`: Demultiplex raw sequencing data.

* `mapping/cellranger_count`: Align fastq files using Cell Ranger count.

* `mapping/cellranger_count_split`: Split 10x Cell Ranger output directory into separate output fields.

Imported workflows from `czbiohub/utilities`:

* `workflows/1_ingestion/cellranger`: Use Cell Ranger to preprocess 10x data.

* `workflows/1_ingestion/cellranger_demux`: Use cellranger demux to demultiplex sequencing BCL output to FASTQ.

* `workflows/1_ingestion/cellranger_mapping`: Use cellranger count to align 10x fastq files to a reference.


## MINOR CHANGES

* Fix `interactive/run_cirrocumulus` script raising `NotImplementedError` caused by using `MutData.var_names_make_unique()` 
on each modality instead of on the whole `MuData` object.

* Fix `transform/normalize_total` and `interactive/run_cirrocumulus` component build missing a hdf5 dependency.

* `interactive/run_cellxgene`: Updated container to ubuntu:focal because it contains python3.6 but cellxgene dropped python3.6 support.

* `mapping/bd_rhapsody_wta`: Set `--parallel` to true by default.

* `mapping/bd_rhapsody_wta`: Translate Bash script into Python.

* `download/sync_test_resources`: Add `--dryrun`, `--quiet`, and `--delete` arguments.

* `convert/from_h5mu_to_seurat`: Use `eddelbuettel/r2u:22.04` docker container in order to speed up builds by downloading precompiled R packages.

* `mapping/cellranger_count`: Use 5Gb for testing (to adhere to github CI runner memory constraints).

* `convert/from_bdrhap_to_h5ad`: change test data to output from `mapping/bd_rhapsody_wta` after reducing the BD Rhapsody test data size.

* Various `config.vsh.yaml`s: Renamed `values:` to `choices:`.

* `download/download_file` and `transfer/publish`: Switch base container from `bash:5.1` to `python:3.10`.

* `mapping/bd_rhapsody_wta`: Make sure procps is installed.

## BUG FIXES

* `mapping/bd_rhapsody_wta`: Use a smaller test dataset to reduce test time and make sure that the Github Action runners do not run out of disk space.

* `download/sync_test_resources`: Disable the use of the Amazon EC2 instance metadata service to make script work on Github Actions runners.

* `convert/from_h5mu_to_seurat`: Fix unit test requiring Seurat by using native R functions to test the Seurat object instead.

* `mapping/cellranger_count` and `bcl_demus/cellranger_mkfastq`: cellranger uses the `--parameter=value` formatting instead of `--parameter value` to set command line arguments.

* `mapping/cellranger_count`: `--nosecondary` is no longer always applied.

* `mapping/bd_rhapsody_wta`: Added workaround for bug in Viash 0.5.12 where triple single quotes are incorrectly escaped (viash-io/viash#139).

## DEPRECATED

* `bcl_demux/cellranger_mkfastq`: Duplicate of `demux/cellranger_mkfastq`.

# openpipeline 0.3.0

* Add `tx_processing` pipeline with following components:
  - `filter_with_counts`
  - `filter_with_scrublet`
  - `filter_with_hvg`
  - `do_filter`
  - `normalize_total`
  - `regress_out`
  - `log1p`
  - `pca`
  - `find_neighbors`
  - `leiden`
  - `umap`

# openpipeline 0.2.0

## NEW FUNCTIONALITY

* Added `from_10x_to_h5ad` and `download_10x_dataset` components.

## MINOR CHANGES
* Workflow `bd_rhapsody_wta`: Minor change to workflow to allow for easy processing of multiple samples with a tsv.

* Component `bd_rhapsody_wta`: Added more parameters, `--parallel` and `--timestamps`.

* Added `pbmc_1k_protein_v3` as a test resource.

* Translate `bd_rhapsody_extracth5ad` from R into Python script.

* `bd_rhapsody_wta`: Remove temporary directory after execution.

* `files/make_params`: Implement unit tests (PR #505).

# openpipeline 0.1.0

* Initial release containing only a `bd_rhapsody_wta` pipeline and corresponding components.
