# ---- libraries ----
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
#library(reticulate)
library(glue)

# ---- load functions ----
source("./r_functions/project_X_to_Y_proteins.R")
source("./r_functions/resultsTabUI.R")
source("./r_functions/inputTabUI.R")
source("./r_functions/examplesTabUI.R")
source("./r_functions/run_python_pipeline.R")
source("./r_functions/plot_figures.R")

# ---- App / paths configuration ----
# APP_DIR <- Sys.getenv("APP_DIR", "/app")
APP_DIR <- Sys.getenv("APP_DIR", tempdir())
APP_DIR = "./"

gsea_from_env <- Sys.getenv("GSEA_PATHWAY_DIR", unset = "")
gsea_pathway_dir <- if (nzchar(gsea_from_env)) {
  normalizePath(gsea_from_env, winslash = "/", mustWork = FALSE)
} else {
  file.path(APP_DIR, "gsea_pathways")
}
dir.create(gsea_pathway_dir, recursive = TRUE, showWarnings = FALSE)
options(gsea_pathway_dir = gsea_pathway_dir)

out_date <- as.character(Sys.Date())
PLOTS_DIR  <- Sys.getenv("PLOTS_DIR", file.path(APP_DIR, "/plots"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", file.path(APP_DIR, "/output"))
PYTHON_OUTPUT_DIR <- Sys.getenv("PYTHON_OUTPUT_DIR", file.path(APP_DIR, "/output/python_output"))

## remove existing output folder
unlink(PYTHON_OUTPUT_DIR, recursive = TRUE, force = TRUE)
unlink(OUTPUT_DIR, recursive = TRUE, force = TRUE)

dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PYTHON_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- UI / Theme (unchanged except small add) ----

## ----bs_theme setup ----

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

## ----inputTabUI ----



## ----placeholder_image_card ----
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

## ----resultsTabUI ----


## ----examplesTabUI ----


## ----fluidPageUI ----
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


# ----SERVER PART----


server <- function(input, output, session){
  results      <- reactiveVal(NULL)
  output_plots <- reactiveVal(character(0))
  out_date_val <- reactiveVal(as.character(Sys.Date()))  # will be updated per run

  options(shiny.maxRequestSize = 500 * 1024^2)

## ----- Upload debug in UI -----
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
    new_dataset <- list(intensity = mat, meta = list(sample = colnames(mat), tissue = NA, diet = NA))
    
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
  

  
  
## ---- processed_data ----
###eventReactive: run the whole pipeline once button pressed
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

    ###ensure each element has intensity matrix and meta
    for (nm in names(data_stes)) {
      if (is.matrix(data_stes[[nm]])) {
        data_stes[[nm]] <- list(intensity = data_stes[[nm]],
                                meta = list(tissue = NA, diet = NA))
      } else {
        if (is.null(data_stes[[nm]]$intensity)) {
          stop("Each dataset must contain an 'intensity' matrix. Problem with: ", nm)
        }
        if (is.null(data_stes[[nm]]$meta)) data_stes[[nm]]$meta <- list(tissue = NA, diet = NA)
      }
    }

#### ---- export meta data ----    
    ## make a new data frame with meta data info
    meta_data = data.frame(
      file_name = character(),
      label = character(),
      tissue = character(),
      diet = character(),
      stringsAsFactors = FALSE
    )
    for (l in names(data_stes)) {
      
      tissue_type <- data_stes[[l]][["meta"]][["tissue"]]
      diet_type <- data_stes[[l]][["meta"]][["diet"]]
      cohort = strsplit(l, "_")[[1]][1]
      meta_data <- rbind(meta_data, data.frame(
        file_name = paste0(cohort,"-",cohort,"-",data_stes[[l]][["meta"]][["sample"]]),
        label = paste0(tissue_type,"-",diet_type),
        tissue = tissue_type,
        diet = diet_type,
        stringsAsFactors = FALSE
      ))
      
    }
    
    meta_data$diet = ifelse(meta_data$diet == "High Fat Diet", "HFD", 
                            ifelse(meta_data$diet %in% c("Standard Diet","AdLib"), "CD", meta_data$diet))
    
    
    write.csv(meta_data, file = paste0(OUTPUT_DIR,"/meta_data", out_date, ".csv"), row.names = F)
    
    
    print("Meta data extraction complete")
    
### ---- normalization ----
    
    # Simple normalization: log1p + column-centering if max > threshold
    max_vals <- sapply(data_stes, function(x) max(x$intensity, na.rm = TRUE))
    datasets_to_normalize <- names(max_vals)[which(max_vals > 0)]

    for(i in names(data_stes)){
      if(i %in% datasets_to_normalize){
        log_data <- log1p(as.matrix(data_stes[[i]]$intensity))
        ## replace inf with NA
        log_data[is.infinite(log_data)] <- NA
        
        med_col <- apply(log_data, 2, function(x) median(x, na.rm = TRUE))
        data_stes[[i]]$normalized <- sweep(log_data, 2, med_col, FUN = "-")
        rownames(data_stes[[i]]$normalized) <- rownames(data_stes[[i]]$intensity)
        colnames(data_stes[[i]]$normalized) <- make.unique(colnames(data_stes[[i]]$normalized))
      } else {
        data_stes[[i]]$normalized <- as.matrix(data_stes[[i]]$intensity)
      }
    }
   
    print("Normalization Complete")
    
## ---- COMBAT based batch correction ----
    
#### ---- integrate all dataset together ----
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
    
    # Make a simple dot plot and save it (robust)
    normalization_figure <- try({
      
      pca_plot(input_data=all_matrices_mat,
               meta_data=meta_data,
               file_name="PCA_after_normalization",
               plot_dir = plot_dir)
      
      ## convert to data frame for ggplot
      all_matrices_mat_df <- as.data.frame(all_matrices_mat)
      all_matrices_mat_df$Protein <- rownames(all_matrices_mat_df)
      all_matrices_mat_long <- pivot_longer(all_matrices_mat_df, cols = -Protein, names_to = "Sample", values_to = "Intensity")
      
      ## add batch info
      all_matrices_mat_long <- all_matrices_mat_long %>%
        mutate(Batch = sapply(strsplit(Sample, "-"), function(x) x[1]))
      
      ## ggplot
      p1 <- ggplot(all_matrices_mat_long, aes(x = Batch, y = Intensity)) +
        geom_boxplot(alpha = 0.3, size = 0.5) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(title = "Normalization Check: Protein Intensities Across Samples", x = "Samples", y = "Intensity")
      
      print(file.path(plot_dir, "normalization_boxplot.png"))
      ggsave(filename = file.path(plot_dir, "normalization_boxplot.png"), plot = p1, width = 8, height = 6, dpi = 150)
      
    }, silent = FALSE)
    if (!inherits(normalization_figure, "try-error")) message("Saved normalization figure")
    
    
    # Keep proteins without all NA columns
    common_prot <- all_matrices_mat[rowSums(is.na(all_matrices_mat)) < ncol(all_matrices_mat), , drop = FALSE]

#### ---- batch data ----
    df_descroption <- data.frame(ID = colnames(all_matrices_mat), sample = 1:ncol(all_matrices_mat))
    batch_vector <- c()
    for (i in 1:length(all_matrices)){
      batch_vector <- c(batch_vector, rep(i, ncol(all_matrices[[i]])))
    }
    df_descroption$batch <- batch_vector
    # write.csv(df_descroption, file = paste0(OUTPUT_DIR,"/all_normalized_description.csv"), row.names = F)
    # batch_data <- read.csv(paste0(OUTPUT_DIR,"/all_normalized_description.csv"), sep = ",", header = TRUE)
    batch_data <- df_descroption

    
#### ---- COMBAT ----
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
    
    # Make a simple PCA plot and save it (robust)
    pca_combat <- try({
      pca_plot(input_data=combat_exprs,
               meta_data=meta_data,
               file_name="PCA_after_COMBAT",
               plot_dir = plot_dir)
      

      
    }, silent = TRUE)
    if (!inherits(pca_combat, "try-error")) message("Saved COMBAT PCA")
    
    
## ---- projection  ----
    # Prepare patterns for projecting back
    patterns_list <- setNames(
      vapply(names(data_stes),
             function(nm) strsplit(nm, "_")[[1]][1],
             FUN.VALUE = character(1)),
      names(data_stes)
    )

    # Projection helper (kept largely same but defensive)

    # Project for each dataset
    for (dataset_name in names(data_stes)) {
      X_count <- as.matrix(data_stes[[dataset_name]][["normalized"]])
      ## add the batch information to the sample names
      colnames(X_count) <- paste0(strsplit(dataset_name, "_")[[1]][1],"-",colnames(X_count))
      
      ## convert infinit to NA
      X_count[is.infinite(X_count)] <- NA
      
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
      message("Projection finished for ", dataset_name)
      }

    print("Projection Complete")
    
    

## ---- imputation ----
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

    # Make a simple PCA plot and save it (robust)
    pca_projection <- try({

      pca_plot(input_data=all_matrices_proj ,
               meta_data=meta_data,
               file_name="PCA_after_projection",
               plot_dir = plot_dir)



    }, silent = TRUE)
    if (!inherits(pca_combat, "try-error")) message("Saved PROJECTION PCA")


    missing_percentage <- rowSums(is.na(all_matrices_proj)) / ncol(all_matrices_proj) * 100
    rows_to_impute <- missing_percentage <= 22
    matrix_to_impute <- as.matrix(all_matrices_proj[rows_to_impute, , drop = FALSE])
    colnames(matrix_to_impute) <- make.unique(colnames(matrix_to_impute))

    # Create MSnSet, impute using QRILC
    msnset_object <- MSnSet(exprs = matrix_to_impute)
    imputed_result <- imputeLCMD::impute.QRILC(exprs(msnset_object))
    imputed_matrix <- imputed_result[[1]]

    print("Imputed Matrix Complete")
    out_date_val(out_date)  # store for download handler

    # Save final outputs (RDS & CSV)
    # out_rds <- file.path(OUTPUT_DIR, paste0("imputed_matrix_", out_date, ".rds"))
    # saveRDS(imputed_matrix, file = out_rds)
    # message("Saved imputed matrix to: ", out_rds)


    out_csv <- file.path(OUTPUT_DIR, paste0("imputed_matrix_", out_date, ".csv"))
    write.csv(imputed_matrix, file = paste0(OUTPUT_DIR,"/imputed_matrix_", out_date, ".csv"), row.names = TRUE)
    message("Saved imputed CSV to: ", out_csv)

    # Make a simple PCA plot and save it (robust)
    pca_ok <- try({
      pr <- prcomp(t(imputed_matrix), center = TRUE, scale. = TRUE)
      pc_df <- as.data.frame(pr$x[, 1:2, drop = FALSE])
      pc_df$sample <- rownames(pc_df)

      pc_df = merge(pc_df, meta_data, by.x = "sample", by.y = "file_name", all.x = TRUE)
      pc_df$batch = strsplit(pc_df$sample, "-") %>% sapply(function(x) x[1])
      write.csv(pc_df , file = paste0(OUTPUT_DIR,"/PCA_of_imputed_matrix_", out_date, ".csv"), row.names = TRUE)
      p1 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = batch)) +
        geom_point() +
        # geom_text(hjust = 1.2, size = 3) +
        ggtitle("PCA_after_imputation (color by batch)")

      p2 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = tissue)) +
        geom_point() +
        # geom_text(hjust = 1.2, size = 3) +
        ggtitle("PCA_after_imputation (color by tissue)")

      p3 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = diet)) +
        geom_point() +
        # geom_text(hjust = 1.2, size = 3) +
        ggtitle("PCA_after_imputation (color by diet)")

      # out_plot <- file.path(plot_dir, "pca_tissue_imputed.png")
      ggsave(filename = file.path(plot_dir, "pca_batch_imputed.png"), plot = p1, width = 8, height = 6, dpi = 150)
      ggsave(filename = file.path(plot_dir, "pca_tissue_imputed.png"), plot = p2, width = 8, height = 6, dpi = 150)
      ggsave(filename = file.path(plot_dir, "pca_diet_imputed.png"), plot = p3, width = 8, height = 6, dpi = 150)

    }, silent = TRUE)
    if (!inherits(pca_ok, "try-error")) message("Saved PCA plot")

    #Run python pipeline if available (non-fatal)
    
    # # Non-blocking attempt to run python (errors are non-fatal)
     try({
      run_python_pipeline(protein_data_path = paste0(OUTPUT_DIR,"/imputed_matrix_", out_date, ".csv"),
                          sample_labels_path = paste0(OUTPUT_DIR,"/meta_data", out_date, ".csv"),
                          log_file = file.path(OUTPUT_DIR, "/browning_pipeline.log")
      )
      }, silent = FALSE)
  
    
## ---- plot machine-learning results ----
    ai_prediction = read.csv(paste0(OUTPUT_DIR,"/python_output/all_sample_scores.csv"))
    query_samples = ai_prediction %>% subset(Batch=="Query")
    
    for(sample in query_samples$X){
      plot_browning_score(sample,ai_prediction,plot_dir)
    }
    
    
## ----Return results----
    list(
      # imputed_matrix = imputed_matrix,
      ai_prediction = NULL,  # placeholder for future
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

## ----Preview table (shows the imputed_matrix head)----
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

## ----Results table----
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

## ----Plots gallery UI----
  
### ---- processed data ----
  
  observeEvent(processed_data(), {
    
    pd <- processed_data()
    req(pd)
    plot_dir <- pd$plot_dir
    req(dir.exists(plot_dir))
    
    # Manually order your 4 plots
    img_files <- c(
      "normalization_boxplot.png",
      "PCA_after_normalization.png",
      "PCA_after_COMBAT.png",
      "PCA_after_projection.png",
      "pca_batch_imputed.png",
      "pca_tissue_imputed.png"
    )
    
    # Only keep files that exist
    imgs_ordered <- img_files[img_files %in% list.files(plot_dir)]
    imgs_ordered <- paste0(plot_dir, "/", imgs_ordered)
    
    # Debug
    message("DEBUG: imgs_ordered = ", paste(imgs_ordered, collapse = ", "))
    
    # Assign each plot
    if(length(imgs_ordered) >= 1) {
      output$plot1 <- renderImage({
        list(src = imgs_ordered[1], contentType = "image/png", width = 300, height = 250)
      }, deleteFile = FALSE)
    }
    
    if(length(imgs_ordered) >= 2) {
      output$plot2 <- renderImage({
        list(src = imgs_ordered[2], contentType = "image/png", width = 300, height = 250)
      }, deleteFile = FALSE)
    }
    
    if(length(imgs_ordered) >= 3) {
      output$plot3 <- renderImage({
        list(src = imgs_ordered[3], contentType = "image/png", width = 300, height = 250)
      }, deleteFile = FALSE)
    }
    
    if(length(imgs_ordered) >= 4) {
      output$plot4 <- renderImage({
        list(src = imgs_ordered[4], contentType = "image/png", width = 300, height = 250)
      }, deleteFile = FALSE)
    }
    
    if(length(imgs_ordered) >= 5) {
      output$plot5 <- renderImage({
        list(src = imgs_ordered[5], contentType = "image/png", width = 300, height = 250)
      }, deleteFile = FALSE)
    }
    
    if(length(imgs_ordered) >= 6) {
      output$plot6 <- renderImage({
        list(src = imgs_ordered[6], contentType = "image/png", width = 300, height = 250)
      }, deleteFile = FALSE)
    }
    
  })
  
### ---- final results ----
  
  observeEvent(processed_data(), {
    pd <- processed_data()
    req(pd)
    plot_dir <- pd$plot_dir
    req(dir.exists(plot_dir))
    
    # ✅ Look for your Browning score result figures
    img_files <- list.files(plot_dir, pattern = "^Browning_score.*\\.png$", full.names = TRUE)
    req(length(img_files) > 0)
    
    # Debug print
    message("DEBUG (final results): Found ", length(img_files), " final result figures.")
    
    # ---- DYNAMIC UI CREATION ----
    output$final_figures_gallery <- renderUI({
      req(img_files)
      n_cols <- 4  # 4 figures per row
      n_imgs <- length(img_files)
      
      # Split figure indices into groups of 4
      img_groups <- split(seq_len(n_imgs), ceiling(seq_len(n_imgs) / n_cols))
      
      # Build rows dynamically
      tagList(
        lapply(img_groups, function(group) {
          fluidRow(
            lapply(group, function(i) {
              column(
                width = 12 / n_cols,  # 4 columns per row → width = 3
                div(
                  class = "card shadow-sm mb-3 p-2",
                  h4(
                    basename(img_files[i]),
                    style = "text-align:center; font-size:13px; font-weight:normal;"
                  ),
                  imageOutput(paste0("final_plot_", i), width = "100%", height = "250px")
                )
              )
            })
          )
        })
      )
    })
    
    # ---- RENDER EACH IMAGE ----
    for (i in seq_along(img_files)) {
      local({
        my_i <- i
        output[[paste0("final_plot_", my_i)]] <- renderImage({
          list(
            src = img_files[my_i],
            contentType = "image/png",
            width = "100%",
            height = "auto"
          )
        }, deleteFile = FALSE)
      })
    }
  })
  

## ----Download processed data----
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
