task MergePicardMet{
  String target_dir
  String met_name
  String uuid
  String output_name
  command{
    gsutil cp -r ${target_dir} ./   
    python pipeline_testing_scripts/merge_picard_mets.py   ${uuid} ${met_name}  ${output_name}
  }

  output{
    File merged_metrics = "${output_name}"
  }
  runtime{
    docker:"gcr.io/broad-dsde-mint-dev/analysis-tools:0.0.2"
    memory:"7.5 GB"
    disks: "local-disk 100 HDD"
  }
}

workflow run_merge{
  String scripts_dir
  String uuid
  Array[String] metrics_names
  scatter(idx in range(length(metrics_names))){
    call MergePicardMet{
      input:
        target_dir = scripts_dir,
        uuid = uuid,
        met_name = metrics_names[idx],
        output_name = metrics_names[idx],
  }
}
}
