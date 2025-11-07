# app.R - corrected
library(shiny)
library(bslib)
library(DT)
library(RColorBrewer)
library(ggplot2)
library(sva)
library(MSnbase)
library(imputeLCMD)
library(dplyr)
library(tidyr)
library(tibble)
library(Biobase)
library(reticulate)
library(glue)

# -------------------------
# App / paths configuration
# -------------------------
# APP_DIR <- Sys.getenv("APP_DIR", "/app")
APP_DIR <- Sys.getenv("APP_DIR", tempdir())

gsea_from_env <- Sys.getenv("GSEA_PATHWAY_DIR", unset = "")
gsea_pathway_dir <- if (nzchar(gsea_from_env)) {
  normalizePath(gsea_from_env, winslash = "/", mustWork = FALSE)
} else {
  file.path(APP_DIR, "gsea_pathways")
}
dir.create(gsea_pathway_dir, recursive = TRUE, showWarnings = FALSE)
options(gsea_pathway_dir = gsea_pathway_dir)

PLOTS_DIR  <- Sys.getenv("PLOTS_DIR", file.path(APP_DIR, "plots"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", file.path(APP_DIR, "output"))
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# UI / Theme (unchanged except small add)
# -------------------------
my_theme <- bs_theme(
  bootswatch = "darkly",
  version = 4,
  primary = "#005f99",
  secondary = "#3caea3",
  base_font = font_google("Roboto"),
  heading_font = font_google("Roboto"),
  code_font = font_google("Source Code Pro")
) %>%
  bs_add_rules(
    "
    .plot-gallery { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 15px; }
    .plot-gallery .shiny-plot-output { margin: auto; }
    .data-table-container { overflow-x: auto; }
  "
  )

inputTabUI <- function() {
  fluidPage(
    fluidRow(
      column(
        width = 4,
        h3("Upload Data"),
        p("Select your input data file below (CSV or RDS format)."),
        fileInput(
          "data_file", "Choose File",
          accept = c(".csv", ".rds"),
          buttonLabel = "Browse...", placeholder = "No file selected"
        ),
        tags$hr(),
        h4("Upload Status"),
        verbatimTextOutput("debug_upload"),
        actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary"),
        br(), br(),
        strong(textOutput("analysis_status"))
      ),
      column(
        width = 8,
        h3("Data Preview"),
        p("After uploading and running analysis, a preview of the data will appear here."),
        DTOutput("preview_table"),
        verbatimTextOutput("non_table_preview")
      )
    )
  )
}

placeholder_image_card <- function(title, file, height = "250px") {
  div(
    class = "card shadow-sm mb-3 p-3",
    h4(title, class = "card-title"),
    img(
      src = file,
      class = "img-fluid rounded",
      style = glue::glue("max-height: {height}; object-fit: cover;")
    )
  )
}

resultsTabUI <- function() {
  fluidPage(
    fluidRow(
      column(
        width = 12,
        h3("Processed Data"),
        downloadButton("download_data", "Download Data", class = "btn-secondary")
      )
    ),
    fluidRow(
      column(
        width = 12,
        h3("Plots"),
        p("Interactive gallery of output plots."),
        div(
          class = "plot-gallery",
          placeholder_image_card("Plot 1", "placeholder_plot1.png"),
          placeholder_image_card("Plot 2", "placeholder_plot2.png"),
          placeholder_image_card("Plot 3", "placeholder_plot3.png")
        )
      )
    )
  )
}

examplesTabUI <- function() {
  fluidPage(
    h2("Examples & Tutorial"),
    p("This section provides example usage, sample data, and instructions."),
    h2("How to Use the App & Data Formatting Guide"),
    p("This page explains exactly what to do, what happens to your data at each step, and how to format your files so the analysis runs smoothly. Keep it open the first time you use the app."),
    tags$details(
      tags$summary(h3("What happens to your data (pipeline)")),
      tags$h4("Setup"),
      tags$ul(
        tags$li("Initializes the environment and prepares input files."),
        tags$li("Captures file metadata (format, sample size, and identifiers) to ensure compatibility and traceability.")
      ),
      tags$h4("Loading and Trimming"),
      tags$ul(
        tags$li("Reads your dataset (CSV or RDS)."),
        tags$li("Removes unnecessary spaces, formatting issues, or mislabeled entries."),
        tags$li("Keeps only relevant samples and protein identifiers to maintain consistency across datasets.")
      ),
      tags$h4("Basic Normalization"),
      tags$ul(
        tags$li("Adjusts raw intensity values to reduce technical variability."),
        tags$li("Places samples on a comparable scale while preserving biological signal.")
      ),
      tags$h4("Batch Correction"),
      tags$ul(
        tags$li("Identifies and corrects batch effects from multiple experimental runs."),
        tags$li("Removes technical biases so differences reflect biological variation rather than experiment-specific effects.")
      ),
      tags$h4("Matrix Projection"),
      tags$ul(
        tags$li("Aligns protein intensities across datasets with different coverage or missing genes."),
        tags$li("Projects all data into a common reference matrix, enabling direct comparisons in the same analytical space.")
      ),
      tags$h4("Imputation"),
      tags$ul(
        tags$li("Estimates missing values using proteomics-specific statistical methods."),
        tags$li("Fills in gaps based on observed patterns, preserving as much biological information as possible.")
      ),
      tags$h4("Output"),
      tags$ul(
        tags$li("Produces a harmonized, analysis-ready data matrix."),
        tags$li("Enables exploration within the app and allows downloading for downstream statistical or machine learning analysis.")
      )
    ),
    tags$details(
      tags$summary(h3("Accepted file types")),
      tags$ul(
        tags$li("CSV (.csv) — preferred for portability."),
        tags$li("RDS (.rds) — preferred when you already have clean R objects with correct column types."),
        tags$li("Max file size: 200 MB")
      )
    ),
    tags$details(
      tags$summary(h3("File Expectations")),
      tags$ul(
        tags$li(strong("Intensity Data Frame:"), " The intensity data frame will have the registered intensity of each gene analyzed. Missing data should preferably be noted as ", code("NA"), ". Row names = gene names; column names = observations."),
        tags$li(strong("Meta Data Frame:"), " This data frame will provide sample metadata. Row names = observations.")
      )
    ),
    tags$details(
      tags$summary(h3("CSV formatting rules")),
      tags$ul(
        tags$li("Header row in the first line with unique column names."),
        tags$li("Delimiter: comma , (other delimiters supported if specified)."),
        tags$li("Encoding: UTF-8."),
        tags$li("Quotes: double quotes \" around fields that contain commas or line breaks."),
        tags$li("Decimal: . (dot)."),
        tags$li("Missing values: prefer empty cells or NA.")
      )
    ),
    tags$details(
      tags$summary(h3("App workflow (tabs & actions)")),
      tags$h4("1) Data Input"),
      tags$ul(
        tags$li("Upload: Choose CSV or RDS."),
        tags$li("Preview: See the first N rows and summary."),
        tags$li("Validate: The app checks types, missingness, duplicates, and date parsing.")
      ),
      tags$h4("2) Analysis"),
      tags$ul(
        tags$li("Select methods/options: Choose analyses relevant to your task."),
        tags$li("Run Analysis: Executes the pipeline and caches results for this session.")
      ),
      tags$h4("3) Results"),
      tags$ul(
        tags$li("Overview: Key metrics and status of the run."),
        tags$li("Tables: Clean dataset preview, summary tables, and model results."),
        tags$li("Charts: Interactive plots; download buttons on each figure.")
      ),
      tags$h4("4) Download"),
      tags$ul(
        tags$li("Processed data: CSV/RDS."),
        tags$li("Figures: PNG/SVG/PDF as available.")
      )
    ),
    tags$details(
      tags$summary(h3("Reproducibility")),
      p("Each run captures: timestamp, input hash, selected options, and package versions.")
    ),
    tags$details(
      tags$summary(h3("Limits & performance")),
      tags$ul(
        tags$li("Large files: streaming preview only; full analysis happens after you click Run Analysis.")
      )
    ),
    tags$details(
      tags$summary(h3("Privacy & security")),
      p("Data is processed in-memory for your session. It is not shared with other users.")
    ),
    h3("Example Input Data"),
    h3("Sample Output"),
    fluidRow(
      column(width = 6, placeholder_image_card("Sample Plot 1", "placeholder_plot1.png")),
      column(width = 6, placeholder_image_card("Sample Plot 2", "placeholder_plot2.png"))
    )
  )
}

ui <- fluidPage(
  theme = my_theme,
  tags$head(tags$style(HTML(".theme-toggle { margin-right: 10px; color: #ffffff; }"))),
  navbarPage(
    title = "AI-BATS",
    id = "main_navbar",
    tabPanel("Data Input", inputTabUI()),
    tabPanel("Results", resultsTabUI()),
    tabPanel("Examples & Tutorial", examplesTabUI()),
    navbarMenu("Settings",
               tabPanel("Appearance",
                        fluidRow(column(width = 12, align = "center", radioButtons("theme_choice", "Theme:", choices = c("Dark", "Light"), inline = TRUE)))
               ))
  )
)

#################
#    SERVER PART
#################

server <- function(input, output, session){
  results      <- reactiveVal(NULL)
  output_plots <- reactiveVal(character(0))
  out_date_val <- reactiveVal(as.character(Sys.Date()))  # will be updated per run

  options(shiny.maxRequestSize = 500 * 1024^2)

  # ----- Upload debug in UI -----
  output$debug_upload <- renderPrint({
    if (is.null(input$data_file)) return("❌ No file uploaded yet")
    df <- input$data_file
    list(
      name     = df$name,
      size     = df$size,
      type     = df$type,
      datapath = df$datapath,
      exists   = file.exists(df$datapath)
    )
  })
  observeEvent(input$data_file, {
    if (is.null(input$data_file)) return()
    message("📥 fileInput changed: ", input$data_file$name)
    message("   → datapath: ", input$data_file$datapath, " (exists: ", file.exists(input$data_file$datapath), ")")
  }, ignoreInit = FALSE)

  # robust loader: if CSV -> wrap into list with single dataset, if RDS and list -> use as-is
  # Robust loader: CSV upload → append to /data/training.rds → return full list
  data_input <- reactive({
    req(input$data_file)
    ext <- tolower(tools::file_ext(input$data_file$name))
    message("📥 Loading file: ", input$data_file$name, " (ext=", ext, ")")
    
    if (ext != "csv") {
      stop("Please upload a CSV file. RDS upload is not supported in this version.")
    }
    
    # Read the uploaded CSV
    df <- read.csv(input$data_file$datapath,
                   stringsAsFactors = FALSE,
                   check.names = FALSE,
                   row.names = 1)
    mat <- as.matrix(df)
    dataset_name <- tools::file_path_sans_ext(basename(input$data_file$name))
    message("✅ CSV loaded: dim=", paste(dim(mat), collapse = " x "), ", dataset_name=", dataset_name)
    
    # Create the structure expected by the pipeline
    new_dataset <- list(intensity = mat, meta = list(tissue = NA, diet = NA))
    
    # Load existing training list or create one
    training_path <- "./data/training.rds"
    if (file.exists(training_path)) {
      message("📂 Loading existing training list from: ", training_path)
      training_list <- readRDS(training_path)
      if (!is.list(training_list)) {
        warning("Existing training.rds is not a list — reinitializing.")
        training_list <- list()
      }
    } else {
      message("🆕 No training.rds found. Creating a new list.")
      training_list <- list()
    }
    
    # Append new dataset (overwrite if name already exists)
    training_list[["Query"]] <- new_dataset
    
    # Save back to RDS
    # saveRDS(training_list, training_path)
    message("💾 Updated training list saved to: ", training_path,
            " (", length(training_list), " datasets total)")
    
    # Return full list to feed into preprocessing
    return(training_list)
  })
  

  
  
  # eventReactive: run the whole pipeline once button pressed
  processed_data <- eventReactive(input$run_analysis, {
    message("▶️ Run Analysis button pressed")

    # create run-specific plot dir
    plot_dir <- tempfile("plots_")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    message(" [run] plot_dir = ", plot_dir)

    data_stes <- data_input()
    # Expect data_stes to be a named list of dataset objects each with $intensity and $meta
    if (!is.list(data_stes) || length(data_stes) == 0) {
      stop("Uploaded data must be a named list or a single intensity matrix (CSV/RDS).")
    }

    # ensure each element has intensity matrix and meta
    for (nm in names(data_stes)) {
      if (is.matrix(data_stes[[nm]])) {
        data_stes[[nm]] <- list(intensity = data_stes[[nm]], meta = list(tissue = NA, diet = NA))
      } else {
        if (is.null(data_stes[[nm]]$intensity)) {
          stop("Each dataset must contain an 'intensity' matrix. Problem with: ", nm)
        }
        if (is.null(data_stes[[nm]]$meta)) data_stes[[nm]]$meta <- list(tissue = NA, diet = NA)
      }
    }

    # Simple normalization: log1p + column-centering if max > threshold
    max_vals <- sapply(data_stes, function(x) max(x$intensity, na.rm = TRUE))
    datasets_to_normalize <- names(max_vals)[which(max_vals > 10000)]

    for(i in names(data_stes)){
      if(i %in% datasets_to_normalize){
        log_data <- log1p(as.matrix(data_stes[[i]]$intensity))
        med_col <- apply(log_data, 2, function(x) median(x, na.rm = TRUE))
        data_stes[[i]]$normalized <- sweep(log_data, 2, med_col, FUN = "-")
        rownames(data_stes[[i]]$normalized) <- rownames(data_stes[[i]]$intensity)
        colnames(data_stes[[i]]$normalized) <- make.unique(colnames(data_stes[[i]]$normalized))
      } else {
        data_stes[[i]]$normalized <- as.matrix(data_stes[[i]]$intensity)
      }
    }

    print("Normalization Complete")

    #### BATCH CORRECTION
    all_genes <- lapply(data_stes, function(x){ rownames(x$intensity) })
    all_genes <- Reduce(union, all_genes)
    ## remove genes name as ""
    all_genes <- all_genes[nchar(all_genes) > 0]

    n <- names(data_stes)
    mat_names <- unlist(lapply(n, function(x) strsplit(x, "_")[[1]][1]))
    all_matrices <- lapply(data_stes, function(x){
      # safe subsetting: ensure rows exist
      m <- matrix(NA, nrow = length(all_genes), ncol = ncol(x$normalized),
                  dimnames = list(all_genes, colnames(x$normalized)))
      common_r <- intersect(rownames(x$normalized), rownames(m))
      m[common_r, ] <- x$normalized[common_r, , drop = FALSE]
      return(m)
    })

    for(i in seq_along(mat_names)){
      colnames(all_matrices[[i]]) <- paste0(mat_names[i],"-",colnames(all_matrices[[i]]))
    }

    # Put all together
    all_matrices_mat <- Reduce(cbind, all_matrices)
    # Keep proteins without all NA columns
    common_prot <- all_matrices_mat[rowSums(is.na(all_matrices_mat)) < ncol(all_matrices_mat), , drop = FALSE]

    df_descroption <- data.frame(ID = colnames(all_matrices_mat), sample = 1:ncol(all_matrices_mat))
    batch_vector <- c()
    for (i in 1:length(all_matrices)){
      batch_vector <- c(batch_vector, rep(i, ncol(all_matrices[[i]])))
    }
    df_descroption$batch <- batch_vector
    write.csv(df_descroption, file = "all_normalized_description.csv", row.names = F)
    batch_data <- read.csv("all_normalized_description.csv", sep = ",", header = TRUE)
    
    batch_names <- c(
      "1" = "Haas2021", 
      "2" = "Haas2022", 
      "3" = "Harney", 
      "4" = "Johanna", 
      "5" = "Kristina ewat", 
      "6" = "Kristina iwat",
      "7" = "Melina",
      "8" = "Oeckl",
      "9" = "Rhabi",
      "10" = "Rosina",
      "11" = "Query"
    )
    
    mat_names <- c("Haas2021", "Haas2022", "Harney", "Johanna", "Kristina", "Kristina", "Melina", "Oeckl", "Rhabi", "Rosina","Query")
    
    all_tissue_types <- c()
    for (l in names(data_stes)) {
      tissue_type <- data_stes[[l]][["meta"]][["tissue"]]
      all_tissue_types <- c(all_tissue_types, tissue_type)
    }
    all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))
    
    # Ensure common_prot is numeric matrix before ComBat
    common_prot <- as.data.frame(common_prot)
    print("Common Proteins Complete")
    
    # Coerce all columns safely
    common_prot[] <- lapply(common_prot, function(col) {
      if (is.character(col) || is.factor(col)) {
        col[col %in% c("NA", "", "NaN")] <- NA
        suppressWarnings(as.numeric(as.character(col)))
      } else {
        as.numeric(col)
      }
    })
    
    # Back to matrix
    common_prot <- as.matrix(common_prot)
    
    # Drop any rows that are entirely NA
    # common_prot <- common_prot[rowSums(is.na(common_prot)) < ncol(common_prot), , drop = FALSE]
    # Drop any rows that contain NA
    common_prot <- common_prot[complete.cases(common_prot), , drop = FALSE]
    
    stopifnot(is.numeric(common_prot))
    
    combat_exprs <- sva::ComBat(
      dat = common_prot,
      batch = batch_data$batch,
      mod = NULL,
      par.prior = TRUE,
      prior.plots = FALSE
    )
    
    print("COMBAT Correction Complete")
    
    #########PROJECTIONS
    # Prepare patterns for projecting back
    patterns_list <- setNames(
      vapply(names(data_stes),
             function(nm) strsplit(nm, "_")[[1]][1],
             FUN.VALUE = character(1)),
      names(data_stes)
    )

    # Projection helper (kept largely same but defensive)
    project_X_to_Y_proteins <- function(X, Y, k = 2, d = 4) {
      X <- as.matrix(X)
      Y <- as.matrix(Y)
      if (is.null(dim(X))) X <- matrix(X, ncol = 1, dimnames = list(names(X), "Sample1"))
      if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1, dimnames = list(names(Y), "Sample1"))

      rn <- rownames(X)
      valid_rows <- !is.na(rn) & nzchar(rn)
      X <- X[valid_rows, , drop = FALSE]
      rownames(X) <- rn[valid_rows]

      common_cols <- intersect(colnames(X), colnames(Y))
      if (length(common_cols) == 0) {
        warning("No common sample columns between X and Y → returning Y as-is")
        return(as.data.frame(Y))
      }
      Xs <- X[, common_cols, drop = FALSE]
      Ys <- Y[, common_cols, drop = FALSE]

      if (ncol(Xs) != ncol(Ys)) stop("Projection aborted: mismatch in sample counts")

      all_missing_proteins <- setdiff(rownames(Xs), rownames(Ys))
      all_missing_proteins <- all_missing_proteins[nchar(all_missing_proteins) > 0]
      all_common_proteins  <- intersect(rownames(Xs), rownames(Ys))
      Y_proj <- Ys

      pr_X <- try(princomp(Xs), silent = TRUE)
      if (inherits(pr_X, "try-error") || is.null(pr_X$scores)) {
        # fallback: use rows as-is
        for (j in all_missing_proteins) {
          avg <- matrix(rowMeans(Xs[j, , drop = FALSE], na.rm = TRUE), nrow = 1)
          rownames(avg) <- j
          colnames(avg) <- colnames(Ys)
          Y_proj <- rbind(Y_proj, avg)
        }
        return(as.data.frame(Y_proj))
      }
      if (is.null(rownames(pr_X$scores))) rownames(pr_X$scores) <- rownames(Xs)
      d_use <- min(d, ncol(pr_X$scores))
      X_dist <- as.matrix(dist(pr_X$scores[, 1:d_use, drop = FALSE], method = "euclidean"))

      for (protein_to_project in all_missing_proteins) {
        if (!protein_to_project %in% rownames(X_dist)) next
        neighbors <- sort(X_dist[protein_to_project, all_common_proteins, drop = FALSE])[1:min(k, length(all_common_proteins))]
        closest_proteins <- names(neighbors)
        if (length(closest_proteins) == 0) next

        batch_vec <- Xs[closest_proteins, , drop = FALSE] - Ys[closest_proteins, , drop = FALSE]
        projected_vec <- matrix(0, nrow = nrow(batch_vec), ncol = ncol(batch_vec))
        for (ii in seq_len(nrow(batch_vec))) {
          projected_vec[ii, ] <- Xs[protein_to_project, ] - batch_vec[ii, ]
        }
        average_proj <- matrix(colMeans(projected_vec, na.rm = TRUE), nrow = 1, ncol = ncol(Ys))
        rownames(average_proj) <- protein_to_project
        colnames(average_proj) <- colnames(Ys)
        Y_proj <- rbind(Y_proj, average_proj)
      }
      return(as.data.frame(Y_proj))
    }

    # Project for each dataset
    for (dataset_name in names(data_stes)) {
      X_count <- as.matrix(data_stes[[dataset_name]][["normalized"]])
      X_count[is.na(X_count)] <- 0
      pattern <- patterns_list[[dataset_name]]
      # select columns from combat_exprs by pattern (case-insensitive)
      cols <- grep(pattern, colnames(combat_exprs), ignore.case = TRUE, value = TRUE)
      if (length(cols) == 0) {
        # fallback: use all columns
        Y_count <- combat_exprs
      } else {
        Y_count <- combat_exprs[, cols, drop = FALSE]
      }
      data_stes[[dataset_name]][["projection"]] <- as.data.frame(project_X_to_Y_proteins(X_count, Y_count))
    }

    print("Projection Complete")

    # Build final projection matrix across datasets
    all_genes <- lapply(data_stes, function(x) rownames(x$projection))
    all_genes <- Reduce(union, all_genes)
    mat_names <- unlist(lapply(names(data_stes), function(x) strsplit(x, "_")[[1]][1]))
    all_matrices <- lapply(seq_along(data_stes), function(i) {
      name_i <- names(data_stes)[i]
      m <- data_stes[[name_i]]$projection
      mat <- matrix(NA, nrow = length(all_genes), ncol = ncol(m), dimnames = list(all_genes, colnames(m)))
      common_r <- intersect(rownames(m), all_genes)
      mat[common_r, ] <- as.matrix(m[common_r, , drop = FALSE])
      colnames(mat) <- paste0(mat_names[i], "-", colnames(mat))
      return(mat)
    })
    all_matrices_proj <- Reduce(cbind, all_matrices)

    missing_percentage <- rowSums(is.na(all_matrices_proj)) / ncol(all_matrices_proj) * 100
    rows_to_impute <- missing_percentage <= 50
    matrix_to_impute <- as.matrix(all_matrices_proj[rows_to_impute, , drop = FALSE])
    colnames(matrix_to_impute) <- make.unique(colnames(matrix_to_impute))

    # Create MSnSet, impute using QRILC
    msnset_object <- MSnSet(exprs = matrix_to_impute)
    imputed_result <- imputeLCMD::impute.QRILC(exprs(msnset_object))
    imputed_matrix <- imputed_result[[1]]

    print("Imputed Matrix Complete")
    out_date <- as.character(Sys.Date())
    out_date_val(out_date)  # store for download handler

    # Save final outputs (RDS & CSV)
    out_rds <- file.path(OUTPUT_DIR, paste0("imputed_matrix_", out_date, ".rds"))
    saveRDS(imputed_matrix, file = out_rds)
    message("Saved imputed matrix to: ", out_rds)
    

    out_csv <- file.path(OUTPUT_DIR, paste0("imputed_matrix_", out_date, ".csv"))
    write.csv(imputed_matrix, file = paste0("imputed_matrix_", out_date, ".csv"), row.names = TRUE)
    message("Saved imputed CSV to: ", out_csv)

    # Make a simple PCA plot and save it (robust)
    pca_ok <- try({
      pr <- prcomp(t(imputed_matrix), center = TRUE, scale. = TRUE)
      pc_df <- as.data.frame(pr$x[, 1:2, drop = FALSE])
      pc_df$sample <- rownames(pc_df)
      p <- ggplot(pc_df, aes(x = PC1, y = PC2, label = sample)) + geom_point() + geom_text(hjust = 1.2, size = 3) + ggtitle("Sample PCA (imputed matrix)")
      out_plot <- file.path(plot_dir, "pca_imputed.png")
      ggsave(filename = out_plot, plot = p, width = 8, height = 6, dpi = 150)
      TRUE
    }, silent = TRUE)
    if (!inherits(pca_ok, "try-error")) message("Saved PCA plot")

    # Run python pipeline if available (non-fatal)
    run_python_pipeline <- function(python_script = "./python_scripts/browning_pipeline.py",
                                    protein_data_path = out_csv,
                                    sample_labels_path = NULL,
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

    # Non-blocking attempt to run python (errors are non-fatal)
    try(run_python_pipeline(protein_data_path = out_csv), silent = TRUE)

    # Return results
    list(
      imputed_matrix = imputed_matrix,
      plot_dir = plot_dir,
      out_date = out_date
    )
  }, ignoreNULL = FALSE)

  # observer that reacts to the button and updates UI using processed_data()
  observeEvent(input$run_analysis, {
    message("⚡ collect processed_data() result and update UI")
    pd <- NULL
    try({
      pd <- processed_data()
    }, silent = FALSE)

    if (is.null(pd)) {
      output$analysis_status <- renderText("❌ processed_data() returned NULL — check console.")
      message("processed_data returned NULL")
      return()
    }
    results(pd)
    output$analysis_status <- renderText("✅ Analysis finished.")
    if (!is.null(pd$plot_dir) && dir.exists(pd$plot_dir)) {
      pngs <- list.files(pd$plot_dir, pattern = "\\.png$", full.names = TRUE)
      output_plots(pngs)
      message("🖼️ Found ", length(pngs), " plot(s) in ", pd$plot_dir)
    } else {
      output_plots(character(0))
    }
  }, ignoreInit = TRUE)

  # Preview table (shows the imputed_matrix head)
  output$preview_table <- renderDT({
    pd <- results()
    req(!is.null(pd))
    if (!is.null(pd$imputed_matrix)) {
      datatable(head(as.data.frame(pd$imputed_matrix), 5), options = list(dom = 't', paging = FALSE))
    } else {
      datatable(data.frame(Note = "No matrix in results"), options = list(dom = 't'))
    }
  })

  output$non_table_preview <- renderPrint({
    obj <- results()
    if (is.null(obj)) {
      "No processed object yet."
    } else {
      str(obj, max.level = 1)
    }
  })

  # Results table
  output$results_table <- renderDT({
    dat <- results(); req(!is.null(dat))
    if (is.data.frame(dat)) {
      datatable(dat, options = list(pageLength = 10, autoWidth = TRUE))
    } else if (is.list(dat) && !is.null(dat$imputed_matrix)) {
      datatable(as.data.frame(dat$imputed_matrix), options = list(pageLength = 10, autoWidth = TRUE))
    } else {
      datatable(data.frame(Note = "No renderable data frame in results"), options = list(dom = 't'))
    }
  })

  # Plots gallery UI
  output$plots_gallery <- renderUI({
    imgs <- output_plots(); req(length(imgs) > 0)
    tagList(lapply(seq_along(imgs), function(i) {
      imgId <- paste0("plot_img_", i)
      output[[imgId]] <- renderImage({
        list(src = imgs[i], contentType = "image/png", width = 300)
      }, deleteFile = FALSE)
      imageOutput(imgId, width = "300px")
    }))
  })

  # Download processed data
  output$download_data <- downloadHandler(
    filename = function() {
      od <- out_date_val()
      paste0("processed_data_", od, ".csv")
    },
    content = function(file) {
      dat <- results(); req(!is.null(dat))
      if (is.list(dat) && !is.null(dat$imputed_matrix)) {
        write.csv(as.data.frame(dat$imputed_matrix), file = file, row.names = TRUE)
      } else if (is.data.frame(dat)) {
        write.csv(dat, file = file, row.names = FALSE)
      } else {
        writeLines("Output is not a data frame.", file)
      }
    }
  )
}

shinyApp(ui = ui, server = server)
