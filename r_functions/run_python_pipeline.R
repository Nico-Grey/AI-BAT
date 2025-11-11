run_python_pipeline <- function(python_script = "./python_scripts/browning_pipeline.py",
                                protein_data_path = paste0(OUTPUT_DIR,"imputed_matrix_", out_date, ".csv"),
                                sample_labels_path = paste0(OUTPUT_DIR,"meta_data", out_date, ".csv"),
                                output_dir = OUTPUT_DIR,
                                log_file = file.path(output_dir, "browning_pipeline.log")) {
  py3 <- Sys.which("python3")
  if (!nzchar(py3)) {
    message("python3 not found on PATH — skipping Python pipeline.")
    return(invisible(NULL))
  }
  python_script_full <- file.path(getwd(), "python_scripts", "browning_pipeline.py")
  if (!file.exists(python_script_full)) {
    message("Python script not found at ", python_script_full, " — skipping Python pipeline.")
    return(invisible(NULL))
  }
  args <- c(python_script_full, "--protein_data_path", protein_data_path)
  if (!is.null(sample_labels_path)) args <- c(args, "--sample_labels_path", sample_labels_path)
  exit_code <- system2(command = py3, args = args, stdout = log_file, stderr = log_file)
  if (exit_code != 0) {
    message("Python pipeline returned non-zero exit code: ", exit_code, ". See log: ", log_file)
    return(invisible(NULL))
  } else {
    message("Python pipeline completed successfully. Log: ", log_file)
    return(invisible(TRUE))
  }
}