nextflow.enable.dsl=2

include { scvi_knn } from params.rootDir + "/target/nextflow/workflows/annotation/scvi_knn/main.nf"
include { scvi_knn_test } from params.rootDir + "/target/nextflow/test_workflows/annotation/scvi_knn_test/main.nf"

workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        obsm_embedding: "X_pca",
        input_obs_batch_label: "sample_id",
        reference_obs_batch_label: "donor_assay",
        reference_var_gene_names: "ensemblid",
        reference_obs_targets: "cell_type",
        leiden_resolution: [1.0, 0.25]
      ],
      [
        id: "no_leiden_resolutions_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        obsm_embedding: "X_pca",
        input_obs_batch_label: "sample_id",
        reference_obs_batch_label: "donor_assay",
        reference_var_gene_names: "ensemblid",
        reference_obs_targets: "cell_type",
        leiden_resolution: []
      ]
    ])
    | map{ state -> [state.id, state] }
    | scvi_knn 
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

    "Output: $output"
    }
    | scvi_knn_test.run(
      fromState: [
        "input": "output"
      ]
    )
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["no_leiden_resolutions_test", "simple_execution_test"]
    }
    }
  