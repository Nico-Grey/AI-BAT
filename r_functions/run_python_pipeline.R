run_python_pipeline <- function(python_script = "./python_scripts/browning_pipeline.py",
                                protein_data_path = paste0("output/","imputed_matrix_", out_date, ".csv"),
                                sample_labels_path = paste0("output/","meta_data", out_date, ".csv"),
                                log_file = file.path("output/", "browning_pipeline.log")) {
  py3 <- Sys.which("python3")

  args <- c(python_script, "--protein_data_path", protein_data_path)
  args <- c(args, "--sample_labels_path", sample_labels_path)

  message("Running Python pipeline starts...")

  system2(command = py3, args = args, stdout = log_file, stderr = log_file)
  ## print the command is running
  message("Command: ", paste(shQuote(py3), paste(shQuote(args), collapse = " ")))

  message("Running Python pipeline ends...")

}