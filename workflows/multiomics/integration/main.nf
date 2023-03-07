nextflow.enable.dsl=2

workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

include { pca } from targetDir + '/dimred/pca/main.nf'
include { find_neighbors } from targetDir + '/neighbors/find_neighbors/main.nf'
include { umap } from targetDir + '/dimred/umap/main.nf'
include { leiden } from targetDir + '/cluster/leiden/main.nf'
include { harmonypy } from targetDir + '/integrate/harmonypy/main.nf'

include { readConfig; helpMessage; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include { setWorkflowArguments; getWorkflowArguments } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/multiomics/integration/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | view { "Input: $it" }
    | run_wf
    | view { "Output: $it" }
}

workflow run_wf {
  take:
  input_ch

  main:
  pca_ch = input_ch
    | preprocessInputs("config": config)
    // split params for downstream components
    | setWorkflowArguments(
      pca: [
        "input": "input", 
        "obsm_output": "obsm_pca",
        "var_input": "var_pca_feature_selection"
      ],
      integration: [
        "obsm_input": "obsm_pca",
        "obs_covariates": "obs_covariates",
        "obsm_output": "obsm_integrated"
      ],
      neighbors: [
        "uns_output": "uns_neighbors",
        "obsp_distances": "obsp_neighbor_distances",
        "obsp_connectivities": "obsp_neighbor_connectivities"
      ],
      clustering: [
        "obsp_connectivities": "obsp_neighbor_connectivities",
        "obs_name": "obs_cluster"
      ],
      umap: [ 
        "uns_neighbors": "uns_neighbors",
        "output": "output",
        "obsm_output": "obsm_umap"
      ]
    )
    | map { tup ->
      assert tup.size() >= 3: "Event should have length 3 or greater."
      def id = tup[0]
      def current_args = tup[1]
      def other_args = tup[2]
      def passthrough = tup.drop(3)
      // If obs_covariates is not set or emtpy, harmony will not be run
      // In this case, the layer that find_neighbour uses should not originate from harmonypy but from pca
      if (!other_args.integration.obs_covariates || 
        (other_args.integration.obs_covariates && other_args.integration.obs_covariates.empty)) {
        other_args.neighbors.put("obsm_input", other_args.pca.obsm_output)
      } else {
        other_args.neighbors.put("obsm_input", other_args.integration.obsm_output)
      }
      [id, current_args, other_args] + passthrough
    }
    | getWorkflowArguments(key: "pca")
    | pca
    | getWorkflowArguments(key: "integration")

  without_harmony_ch = pca_ch
    | filter{!it[1].obs_covariates || (it[1].obs_covariates && it[1].obs_covariates.empty)}

  with_harmony_ch = pca_ch
    | filter{it[1].obs_covariates && !it[1].obs_covariates.empty}
    | harmonypy

  output_ch = without_harmony_ch.concat(with_harmony_ch)
    | getWorkflowArguments(key: "neighbors")
    | find_neighbors

    | getWorkflowArguments(key: "clustering")
    | leiden

    | getWorkflowArguments(key: "umap")
    | umap.run(
      auto: [ publish: true ]
    )

    // remove splitArgs
    | map { tup ->
      tup.take(2) + tup.drop(3)
    }

  emit:
  output_ch
}

workflow test_wf {
  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
      [
        id: "foo",
        input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
        layer: "",
        obs_covariates: "sample_id"
      ],
      [
        id: "foo_without_harmony",
        input: params.resources_test + "/concat_test_data/concatenated_brain_filtered_feature_bc_matrix_subset.h5mu",
        layer: "",
        obs_covariates: []
      ],
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
    | view { "Input: $it" }
    | run_wf
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, file]"
      assert output[1].toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output_list[1]}"
      "Output: $output"
    }
    | toList()
    | map { output_list ->
      println "output_list: $output_list"
      assert output_list.size() == 2 : "output channel should contain two events"
      assert (output_list.collect({it[0]}) as Set).equals(["foo_without_harmony", "foo"] as Set): "Output ID should be same as input ID"
    }
    //| check_format(args: {""}) // todo: check whether output h5mu has the right slots defined
}