workflow run_wf {
  take:
  input_ch

  main:
  reference_ch = input_ch
    // Make sure there is not conflict between the output from this workflow
    // And the output from any of the components
    | map {id, state ->
      def new_state = state + ["workflow_output": state.output]
      [id, new_state]
    }
    // // download the reference h5ad file
    // | download_file.run(
    //   fromState: { id, state ->
    //     [
    //       "input": state.reference_url,
    //       "output": "reference.h5ad",
    //       "verbose": "true",
    //     ]
    //   },
    //   toState: [
    //     "input": "output",
    //   ]
    // )
    // // convert the reference h5ad file to h5mu
    // | from_h5ad_to_h5mu.run(
    //     fromState: { id, state ->
    //     [
    //       "input": state.input,
    //       "modality": "rna",
    //     ]
    //   },
    //   toState: [
    //     "input": "output",
    //   ]
    // )
    | split_samples.run(
        fromState: { id, state ->
        [
          "input": state.reference,
          "modality": "rna",
          "obs_feature": state.obs_reference_batch
        ]
      },
      toState:  ["output": "output", "output_files": "output_files"]
    )

    | flatMap {id, state ->
        def outputDir = state.output
        def files = readCsv(state.output_files.toUriString())
        files.collect{ dat ->
          // def new_id = id + "_" + dat.name
          def new_id = id // it's okay because the channel will get split up anyways
          def new_data = outputDir.resolve(dat.filename)
          [ new_id, state + ["reference_input": new_data]]
        }
        }
    
    // Remove arguments from split samples from state

    | niceView()

    | map {id, state -> 
      def keysToRemove = ["output_files"]
      def newState = state.findAll{it.key !in keysToRemove}
      [id, newState]
    }

    | niceView()

    | view {"After splitting samples: $it"}

    | process_samples_workflow.run(
      fromState: {id, state ->
        def newState = [
          "input": state.reference_input, 
          "id": id,]
      },
      toState: ["processed_reference": "output"]
    )

    | niceView()

    | view {"After splitting samples: $it"}

  output_ch = reference_ch

  emit:
  output_ch
}