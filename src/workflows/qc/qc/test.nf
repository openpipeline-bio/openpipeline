nextflow.enable.dsl=2

include { qc } from params.rootDir + "/target/nextflow/workflows/qc/qc/qc/main.nf"


workflow test_wf {
  // allow changing the resources_test dir
  resources_test = file("${params.rootDir}/resources_test")

  output_ch = 
    Channel.fromList([
      [
        id: "mouse_test",
        input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      ],
      [
        id: "human_test",
        input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      ]
    ])
    | map{ state -> [state.id, state] }
    | qc
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test")

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"

      "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["human_test", "mouse_test"]
    }
  
}