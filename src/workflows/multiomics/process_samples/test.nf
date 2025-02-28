nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow"

include { remove_modality }  from targetDir + '/filter/remove_modality/main.nf'
include { move_layer } from targetDir + '/transform/move_layer/main.nf' 
include { process_samples } from targetDir + "/workflows/multiomics/process_samples/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "mouse",
      input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      rna_min_counts: 2,
      output: "test.h5mu",
    ],
    [
      id: "human",
      input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
      publish_dir: "foo/",
      rna_min_counts: 2,
      output: "test.h5mu",

    ]
  ])
  | map{ state -> [state.id, state] }
  | process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith("test.h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
  
}

workflow test_wf2 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
        id: "pbmc",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        var_name_mitochondrial_genes: 'mitochondrial',
        var_name_ribosomal_genes: 'ribosomal',
        rna_min_counts: 2,
        prot_min_counts: 3,
        add_id_to_obs: true,
        add_id_make_observation_keys_unique: true,
        add_id_obs_output: "sample_id"
    ],
    [
      id: "pbmc_with_more_args",
      input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
      rna_min_counts: 2,
      rna_max_counts: 1000000,
      rna_min_genes_per_cell: 1,
      rna_max_genes_per_cell: 1000000,
      rna_min_cells_per_gene: 1,
      rna_min_fraction_mito: 0.0,
      rna_max_fraction_mito: 1.0,
      prot_min_counts: 3,
      prot_max_counts: 1000000,
      prot_min_proteins_per_cell: 1,
      prot_max_proteins_per_cell: 1000000,
      prot_min_cells_per_protein: 1,
      var_name_mitochondrial_genes: 'mitochondrial',
      obs_name_mitochondrial_fraction: 'fraction_mitochondrial',
      var_name_ribosomal_genes: 'ribosomal',
      obs_name_ribosomal_fraction: 'fraction_ribosomal',
      add_id_to_obs: true,
      add_id_make_observation_keys_unique: true,
      add_id_obs_output: "sample_id"
    ],
  ])
  | map{ state -> [state.id, state] }
  | process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    // The result of this pipeline is always 1 merged sample, regardless of the number of input samples. 
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
}

workflow test_wf3 {

  resources_test = file(params.resources_test)

  input_ch = Channel.fromList([
      [
        id: "mouse",
        input: resources_test.resolve("concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
        publish_dir: "foo/",
        rna_min_counts: 2,
        var_qc_metrics: "highly_variable",
        highly_variable_features_var_output: "highly_variable",
      ],
      [
        id: "human",
        input: resources_test.resolve("concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"),
        publish_dir: "foo/",
        rna_min_counts: 2,
        var_qc_metrics: "highly_variable",
        highly_variable_features_var_output: "highly_variable",
      ]
    ])
    | map{ state -> [state.id, state] }

    human_ch = input_ch
      | filter{it[0] == "human"}
      | remove_modality.run(
        fromState: { id, state ->
          [
            "input": state.input,
            "modality": "atac"
          ]
        },
        toState: ["input": "output"]
      )

    mouse_ch = input_ch
      | filter{it[0] == "mouse"}
      | remove_modality.run(
        fromState: { id, state ->
          [
            "input": state.input,
            "modality": "rna"
          ]
        },
        toState: ["input": "output"]
      )

    output_ch_test_2 = human_ch.concat(mouse_ch)
      | process_samples
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
      }
}

workflow test_wf4 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
        id: "10x_5k_lung_crispr",
        input: resources_test.resolve("10x_5k_lung_crispr/SC3_v3_NextGem_DI_CRISPR_A549_5K.h5mu"),
        rna_min_counts: 2,
        gdo_min_counts: 10,
        gdo_max_counts: 10000,
        add_id_to_obs: true,
        add_id_make_observation_keys_unique: true,
        add_id_obs_output: "sample_id"
    ]
  ])
  | map{ state -> [state.id, state] }
  | process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    // The result of this pipeline is always 1 merged sample, regardless of the number of input samples. 
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
}

workflow test_wf5 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "5k_human_antiCMV_T_TBNK_select_layer",
      input: resources_test.resolve("10x_5k_lung_crispr/SC3_v3_NextGem_DI_CRISPR_A549_5K.h5mu"),
      rna_min_counts: 2,
      rna_max_counts: 1000000,
      rna_min_genes_per_cell: 1,
      rna_max_genes_per_cell: 1000000,
      rna_min_cells_per_gene: 1,
      rna_min_fraction_mito: 0.0,
      rna_max_fraction_mito: 1.0,
      rna_min_fraction_ribo: 0.0,
      rna_max_fraction_ribo: 1.0,
      var_name_mitochondrial_genes: 'mitochondrial',
      obs_name_mitochondrial_fraction: 'fraction_mitochondrial',
      var_name_ribosomal_genes: 'ribosomal',
      obs_name_ribosomal_fraction: 'fraction_ribosomal',
      add_id_to_obs: true,
      add_id_make_observation_keys_unique: true,
      add_id_obs_output: "sample_id",
      rna_layer: "test_layer"
    ],
  ])
  | map{ state -> [state.id, state] }
  | move_layer.run(
    fromState: { id, state ->
      [
        "input": state.input,
        "output_layer": "test_layer",
        "modality": "rna"
      ]
    },
    toState: ["input": "output"]
  )
  | process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    // The result of this pipeline is always 1 merged sample, regardless of the number of input samples. 
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
}

workflow test_wf6 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
        id: "pbmc_clr_axis_0",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        clr_axis: 1
    ],
    [
        id: "pbmc_clr_axis_1",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        clr_axis: 1
    ]
  ])
  | map{ state -> [state.id, state] }
  | process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
  
}


// The following test is supposed to fail. It is used to test the error handling of the pipeline.
// However, there is not way to catch the error originating from the workflow

// workflow test_wf4 {

//   helpMessage(config)

//   // allow changing the resources_test dir
//   params.resources_test = params.rootDir + "/resources_test"
  

//   testParams = [
//     param_list: [
//         [
//           id: "pbmc",
//           input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
//           output: "foo.h5mu"
//         ],
//         [
//           id: "pbmc_with_more_args",
//           input: params.resources_test + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
//           output: "foo2.h5mu"
//         ],
//       ],

//     ]

//     output_ch = channelFromParams(testParams, config)
//       | run_wf
// }

workflow test_wf7 {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
    [
      id: "5k_human_antiCMV_T_TBNK_select_layer",
      input: resources_test.resolve("10x_5k_lung_crispr/SC3_v3_NextGem_DI_CRISPR_A549_5K.h5mu"),
      rna_min_counts: 2,
      rna_max_counts: 1000000,
      rna_min_genes_per_cell: 1,
      rna_max_genes_per_cell: 1000000,
      rna_min_cells_per_gene: 1,
      rna_min_fraction_mito: 0.0,
      rna_max_fraction_mito: 1.0,
      rna_min_fraction_ribo: 0.0,
      rna_max_fraction_ribo: 1.0,
      rna_enable_scaling: true,
      var_name_mitochondrial_genes: 'mitochondrial',
      obs_name_mitochondrial_fraction: 'fraction_mitochondrial',
      var_name_ribosomal_genes: 'ribosomal',
      obs_name_ribosomal_fraction: 'fraction_ribosomal',
      add_id_to_obs: true,
      add_id_make_observation_keys_unique: true,
      add_id_obs_output: "sample_id",
      rna_layer: "test_layer",

    ],
  ])
  | map{ state -> [state.id, state] }
  | move_layer.run(
    fromState: { id, state ->
      [
        "input": state.input,
        "output_layer": "test_layer",
        "modality": "rna"
      ]
    },
    toState: ["input": "output"]
  )
  | process_samples
  | view { output ->
    assert output.size() == 2 : "outputs should contain two elements; [id, file]"
    assert output[1].output.toString().endsWith(".h5mu") : "Output file should be a h5mu file. Found: ${output[1].output}"
    "Output: $output"
  }
  | toSortedList()
  | map { output_list ->
    // The result of this pipeline is always 1 merged sample, regardless of the number of input samples. 
    assert output_list.size() == 1 : "output channel should contain one event"
    assert output_list[0][0] == "merged" : "Output ID should be 'merged'"
  }
}
