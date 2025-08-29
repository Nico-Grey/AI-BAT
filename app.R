# app.R - self-contained (no general_functions.R)
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

# -------------------------
# App / paths configuration
# -------------------------
APP_DIR <- Sys.getenv("APP_DIR", "/app")

# GSEA dir still optional, in case you later want to add it
gsea_from_env <- Sys.getenv("GSEA_PATHWAY_DIR", unset = "")
gsea_pathway_dir <- if (nzchar(gsea_from_env)) {
  normalizePath(gsea_from_env, winslash = "/", mustWork = FALSE)
} else {
  file.path(APP_DIR, "gsea_pathways")
}
dir.create(gsea_pathway_dir, recursive = TRUE, showWarnings = FALSE)
options(gsea_pathway_dir = gsea_pathway_dir)

# Output dirs
PLOTS_DIR  <- Sys.getenv("PLOTS_DIR", file.path(APP_DIR, "plots"))
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR", file.path(APP_DIR, "output"))
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# UI / Theme
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
    /* Gallery for responsive plot layout */
    .plot-gallery {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
      gap: 20px;
      margin-top: 15px;
    }
    .plot-gallery .shiny-plot-output {
      margin: auto;
    }
    /* Ensure table is responsive */
    .data-table-container {
      overflow-x: auto;
    }
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
        # --- Debug block: shows exactly what Shiny received ---
        tags$hr(),
        h4("Upload Status"),
        verbatimTextOutput("debug_upload"),
        tags$hr(),
        selectInput("plot_color", "Choose Point Color:",
                    choices = c("blue", "red", "green", "purple", "black"),
                    selected = "blue"),
        sliderInput("plot_size", "Point Size:", min = 1, max = 10, value = 3),
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

# -------------------------
# Helpers
# -------------------------

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

# -------------------------
# Results Tab (with image placeholders)
# -------------------------

resultsTabUI <- function() {
  fluidPage(
    fluidRow(
      column(
        width = 12,
        h3("Processed Data"),
        downloadButton(
          "download_data",
          "Download Data",
          class = "btn-secondary"
        )
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

# -------------------------
# Examples Tab (with image placeholders)
# -------------------------

examplesTabUI <- function() {
  fluidPage(
    h2("Examples & Tutorial"),
    p("This section provides example usage, sample data, and instructions."),
    
    # Main guide header
    h2("How to Use the App & Data Formatting Guide"),
    p("This page explains exactly what to do, what happens to your data at each step, and how to format your files so the analysis runs smoothly. Keep it open the first time you use the app."),
    
    # Pipeline (collapsible)
    tags$details(
      tags$summary(h3("What happens to your data (pipeline)")),
      tags$h4("Ingest"),
      tags$ul(
        tags$li("Reads CSV (via readr) or RDS (a saved R object, e.g., data.frame/tibble)."),
        tags$li("Captures file metadata (name, size).")
      ),
      
      tags$h4("Validation"),
      tags$ul(
        tags$li("Checks encoding, header row, duplicate column names, and row count."),
        tags$li("Verifies required/selected columns exist and have compatible types."),
        tags$li("Normalizes missing values (e.g., \"\", NA, NaN, NULL).")
      ),
      
      tags$h4("Cleaning"),
      tags$ul(
        tags$li("Trims whitespace, standardizes date/time formats, coerces types (numeric/factor/date), optional de-duplication."),
        tags$li("Flags out-of-range values and uncommon categories.")
      ),
      
      tags$h4("Transformations (only if selected)"),
      tags$ul(
        tags$li("Filtering/subsetting, aggregations, scaling/normalization, encoding of categorical variables, derived fields.")
      ),
      
      tags$h4("Analysis"),
      tags$ul(
        tags$li("Runs the selected methods and diagnostics."),
        tags$li("Computes summary stats and model outputs (as applicable).")
      ),
      
      tags$h4("Results"),
      tags$ul(
        tags$li("Renders interactive tables/plots and downloadable artifacts."),
        tags$li("Data retention: By default, data is kept only in your active session and cleared when you close/refresh the app.")
      )
    ),
    
    # Accepted file types (collapsible)
    tags$details(
      tags$summary(h3("Accepted file types")),
      tags$ul(
        tags$li("CSV (.csv) — preferred for portability."),
        tags$li("RDS (.rds) — preferred when you already have clean R objects with correct column types."),
        tags$li("Max file size: 200 MB")
      )
    ),
    
    # File expectations
    tags$details(
      tags$summary(h3("File Expectations")),
      tags$ul(
        tags$li(
          strong("Intensity Data Frame:"), 
          " The intensity data frame will have the registered intensity of each gene analyzed. ",
          "Missing data should preferably be noted as ", code("NA"), 
          ", and not NaN or 0. The row names of the data frame will be the exact name of the gene analyzed, ",
          "while the column names will be the observations."
        ),
        tags$li(
          strong("Meta Data Frame:"), 
          " This data frame will provide all the additional information for each sample. ",
          "The row names will be the observations (opposite of the Intensity Data Frame). ",
          "The column names need to be written exactly as follows:",
          tags$ul(
            tags$li(code("sample")),
            tags$li(code("treatment")),
            tags$li(code("tissue")),
            tags$li(code("age")),
            tags$li(code("gender")),
            tags$li(code("temperature")),
            tags$li(code("genotype")),
            tags$li(code("mouse"))
          ),
          
          h3("Sample Meta Data Frame"),
          fluidRow(
            column(
              width = 6,
              placeholder_image_card("Meta Data", "example.tab.png")
              
            ),
        )
      )
      )
    ),
    
    # CSV formatting rules
    tags$details(
      tags$summary(h3("CSV formatting rules")),
      tags$ul(
        tags$li("Header row in the first line with unique column names."),
        tags$li("Delimiter: comma , (other delimiters like ; or tab are supported if you specify them in the upload options)."),
        tags$li("Encoding: UTF-8."),
        tags$li("Quotes: double quotes \" around fields that contain commas or line breaks."),
        tags$li("Decimal: . (dot). If you use , (comma), set that option during upload."),
        tags$li("Missing values: prefer empty cells or NA. The app also treats NaN, NULL, and . as missing.")
      )
    ),
    
    # App workflow
    tags$details(
      tags$summary(h3("App workflow (tabs & actions)")),
      tags$h4("1) Data Input"),
      tags$ul(
        tags$li("Upload: Choose CSV or RDS."),
        tags$li("Optional: set delimiter, decimal mark, and NA strings (CSV)."),
        tags$li("Preview: See the first N rows and summary."),
        tags$li("Validate: The app checks types, missingness, duplicates, and date parsing."),
        tags$li("Map columns: Use dropdowns to assign roles (ID, Date, Outcome, Group, Features). Nothing is permanently changed in your file."),
        tags$li("Fix issues: If validation flags problems, you’ll get clear messages and suggested fixes.")
      ),
      
      tags$h4("2) Analysis"),
      tags$ul(
        tags$li("Select methods/options: Choose analyses relevant to your task (summary stats, comparisons, modeling, time series, etc.)."),
        tags$li("Configure settings: Filters, transformations, handling of missing data, resampling, metrics."),
        tags$li("Run Analysis: Executes the pipeline and caches results for this session.")
      ),
      
      tags$h4("3) Results"),
      tags$ul(
        tags$li("Overview: Key metrics and status of the run."),
        tags$li("Tables: Clean dataset preview, summary tables, and model results."),
        tags$li("Charts: Interactive plots; hover to see values; download buttons on each figure."),
        tags$li("Diagnostics: Warnings, assumptions checks, and logs.")
      ),
      
      tags$h4("4) Download"),
      tags$ul(
        tags$li("Processed data: The cleaned/filtered dataset used in the analysis (CSV/RDS)."),
        tags$li("Figures: PNG/SVG/PDF as available."),
        tags$li("Run Report: YAML/JSON describing input file, mappings, options, package versions, random seed.")
      )
    ),
    
    # Reproducibility
    tags$details(
      tags$summary(h3("Reproducibility")),
      p("Each run captures: timestamp, input hash, selected options, and package versions. You can re-run the same configuration by reloading the Run Report (when supported) or by reapplying the settings.")
    ),
    
    # Limits & performance
    tags$details(
      tags$summary(h3("Limits & performance")),
      tags$ul(
        tags$li("Large files: streaming preview only; full analysis happens after you click Run Analysis."),
        tags$li("Memory/timeouts: very wide data or heavy models may take time or fail — reduce columns/rows or simplify options."),
        tags$li("Non-standard CSVs (exotic encodings, mixed delimiters) may need pre-cleaning.")
      )
    ),
    
    # Privacy & security
    tags$details(
      tags$summary(h3("Privacy & security")),
      p("Data is processed in-memory for your session. It is not shared with other users."),
      p("No files are permanently stored unless your deployment explicitly enables it. Replace this line if your server persists data.")
    ),
    
    # Example Input Data & Sample Output (kept from original)
    h3("Example Input Data"),
    p("Replace this with a description or snippet of your example input data format."),
    
    h3("Sample Output"),
    fluidRow(
      column(
        width = 6,
        placeholder_image_card("Sample Plot 1", "placeholder_plot1.png")
      ),
      column(
        width = 6,
        placeholder_image_card("Sample Plot 2", "placeholder_plot2.png")
      )
    )
  )
}



ui <- fluidPage(
  theme = my_theme,
  tags$head(
    tags$style(HTML(".theme-toggle { margin-right: 10px; color: #ffffff; }"))
  ),
  navbarPage(
    title = "AI-BATS",
    id = "main_navbar",
    tabPanel("Data Input", inputTabUI()),
    tabPanel("Results", resultsTabUI()),
    tabPanel("Examples & Tutorial", examplesTabUI()),
    navbarMenu("Settings",
               tabPanel("Appearance",
                        fluidRow(
                          column(
                            width = 12,
                            align = "center",
                            radioButtons("theme_choice", "Theme:", choices = c("Dark", "Light"), inline = TRUE)
                          )
                        )
               )
    )
  )
)


# -------------------------
# Server
# -------------------------
server <- function(input, output, session){

  # reactive holders
  results      <- reactiveVal(NULL)
  output_plots <- reactiveVal(character(0))
  
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
  
  # also log to console when upload changes
  observeEvent(input$data_file, {
    if (is.null(input$data_file)) return()
    message("📥 fileInput changed: ", input$data_file$name)
    message("   → datapath: ", input$data_file$datapath, " (exists: ", file.exists(input$data_file$datapath), ")")
  }, ignoreInit = FALSE)
  
  data_input <- reactive({
    req(input$data_file)
    ext <- tolower(tools::file_ext(input$data_file$name))
    message("📥 fileInput changed: ", input$data_file$name, " (ext=", ext, ")")
    
    if (ext == "rds") {
      dat <- readRDS(input$data_file$datapath)
      message("✅ RDS loaded: class=", paste(class(dat), collapse = ", "), ", length=", ifelse(is.list(dat), length(dat), NA))
      dat
    } else {
      dat <- read.csv(input$data_file$datapath, stringsAsFactors = FALSE)
      message("✅ CSV loaded: dim=", paste(dim(dat), collapse=" x "))
      dat
    }
  })
  
  # --- Run pipeline only when button is pressed ---
  processed_data <- eventReactive(input$run_analysis, {
    message("▶️ Run Analysis button pressed")
    
    # create run-specific plot dir
    plot_dir <- tempfile("plots_")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    message(" [run] plot_dir = ", plot_dir)
    
    data_stes <- data_input()
    validate(
      need(is.list(data_stes), "Uploaded .rds must be a named list of datasets."),
      need(length(data_stes) > 0, "The list of datasets is empty.")
    )
    
    # Check normalization
    max_vals <- c()
    for(i in 1:length(data_stes)){
      max_vals <- c(max_vals, max(data_stes[[i]]$intensity, na.rm = TRUE))
    }
    
    
    datasets_to_normalize <- names(data_stes)[which(max_vals > 10000)]
    
    # Normalize each data set individually
    for(i in names(data_stes)){
      if(i %in% datasets_to_normalize){
        # Calculate the median per column
        log_data <- log1p(data_stes[[i]]$intensity)
        med_col <- apply(log_data, 2, function(x) median(x, na.rm = T))
        
        data_stes[[i]]$normalized <- (log_data - matrix(rep(med_col, nrow(log_data)), nrow = nrow(log_data), byrow = T))
        
        rownames(data_stes[[i]]$normalized) <- rownames(data_stes[[i]]$intensity)
        
      } else {
        # If no normalization is needed, just copy the original data
        data_stes[[i]]$normalized <- data_stes[[i]]$intensity
      }
    }
    print("Normalization Complete")
    
    
    ##Batch Correction and Data Plotting
    all_genes <- lapply(data_stes, function(x){
      return(rownames(x$intensity))
    })
    
    all_genes <- Reduce(union, all_genes)
    
    n <- names(data_stes)
    mat_names <- unlist(lapply(n, function(x) strsplit(x, "_")[[1]][1]))
    all_matrices <- lapply(data_stes, function(x) {
      rn <- rownames(x$normalized)
      
      # Use match to preserve full all_genes length (missing = NA)
      m <- x$normalized[match(all_genes, rn), , drop = FALSE]
      rownames(m) <- all_genes
      
      # Convert factors/characters to numeric safely
      m <- as.data.frame(m)
      m[] <- lapply(m, function(col) {
        if (is.character(col) || is.factor(col)) {
          col[col %in% c("NA", "", "NaN")] <- NA  # unify NA-like strings
          as.numeric(col)
        } else {
          col
        }
      })
      
      return(as.matrix(m))
    })
    
    
    for(i in 1:length(mat_names)){
      colnames(all_matrices[[i]]) <- paste0(mat_names[i],"-",colnames(all_matrices[[i]]))
    }
    
    all_matrices_mat <- Reduce(cbind, all_matrices)
    
    common_prot <- all_matrices_mat[complete.cases(all_matrices_mat), ]
    
    print("Common Protein Created")
    
    
    df_descroption <- data.frame(ID = colnames(all_matrices_mat), sample = 1:ncol(all_matrices_mat))
    batch_vector <- c()
    for (i in 1:length(all_matrices)){
      batch_vector <- c(batch_vector, rep(i, ncol(all_matrices[[i]])))
    }
    df_descroption$batch <- batch_vector
    write.csv(df_descroption, file = file.path(getwd(), "all_normalized_description.csv"), row.names = FALSE)
    batch_data <- read.csv(file.path(getwd(), "all_normalized_description.csv"), sep = ",", header = TRUE)
    
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
      "10" = "Rosina"
    )
    
    mat_names <- c("Haas2021", "Haas2022", "Harney", "Johanna", "Kristina", "Kristina", "Melina", "Oeckl", "Rhabi", "Rosina")
    
    
    print("Batch and Matrices created")
    
    # Loop through each dataset in data_stes
    for (dataset_name in names(data_stes)) {
      # Check if 'meta' exists and contains the column 'protein_type'
      if ("meta" %in% names(data_stes[[dataset_name]]) &&
          "protein_type" %in% colnames(data_stes[[dataset_name]]$meta)) {
        # Rename 'protein_type' to 'tissue'
        colnames(data_stes[[dataset_name]]$meta)[
          colnames(data_stes[[dataset_name]]$meta) == "protein_type"
        ] <- "tissue"
      }
    }
    
    
    all_tissue_types <- c()
    for (l in names(data_stes)) {
      tissue_type <- data_stes[[l]][["meta"]][["tissue"]]
      all_tissue_types <- c(all_tissue_types, tissue_type)
    }
    all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))
    
    print("Tissues Completed")
    
    
    # Perform PCA - Before Batch Correction
    pca_before <- prcomp(t(common_prot), scale.=TRUE)
    pca_data_before <- data.frame(
      PC1 = pca_before$x[, 1], 
      PC2 = pca_before$x[, 2], 
      Batch = factor(batch_data$batch, levels = names(batch_names), labels = mat_names)
    )
    palette <- brewer.pal(n = length(unique(batch_data$batch)), name = "Set3")
    # 1) Build the plot and store it in a variable
    p_before <- ggplot(pca_data_before, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = palette) +
      ggtitle("PCA Plot Before Batch Correction") +
      theme_minimal() +
      theme(
        plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)
      ) +
      labs(color = "Batch")
    
    # 2) Save it into your mounted folder
    ggsave(
      filename = file.path(plots_dir, "pca_before_batch_correction.png"),
      plot     = p_before,
      width    = 8,      # inches → 800 px at dpi=100
      height   = 6,      # inches → 600 px at dpi=100
      units    = "in",
      dpi      = 100
    )
    
    
    
    combat_exprs <- ComBat(dat=common_prot, batch=batch_data$batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
    
    print("COMBAT Completed")
    # Perform PCA - After Batch Correction
    pca_after <- prcomp(t(combat_exprs), scale.=TRUE)
    pca_data_after <- data.frame(
      PC1 = pca_after$x[, 1], 
      PC2 = pca_after$x[, 2], 
      Batch = factor(batch_data$batch, levels = names(batch_names), labels = batch_names)
    )
    pca_after <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +  
      scale_color_manual(values = palette) +
      ggtitle("PCA Plot After Batch Correction") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      ) +
      labs(color = "Batch")
    
    ggsave(
      filename = file.path(plots_dir, "pca_after_batch_correction.png"),
      plot     = pca_after,
      width    = 8,      # inches → 800 px at dpi=100
      height   = 6,      # inches → 600 px at dpi=100
      units    = "in",
      dpi      = 100
    )
    print("Batch Correction Complete")
    
    ##Projections
    # auto-generate a named character vector of “prefixes”
    patterns_list <- setNames(
      vapply(names(data_stes),
             function(nm) strsplit(nm, "_")[[1]][1],
             FUN.VALUE = character(1)),
      names(data_stes)
    )
    
    
    
    #Function for projection
    project_X_to_Y_proteins <- function(X, Y, k = 2, d = 4){
      
      # X: Uncorrected matrix with proteins in the rows and samples in the columns
      # Y: Corrected matrix with proteins in the rows and samples in the columns
      # number of proteins in X should be larger than number of proteins in Y
      # k: the number of nearest neighbors to be considered
      # d: Number of PCs to use
      #Select all the proteins that need to be projected
      
      all_missing_proteins <- setdiff(rownames(X), rownames(Y))
      name_lengths <- nchar(all_missing_proteins)
      all_missing_proteins <- all_missing_proteins[name_lengths > 0]
      all_common_proteins <- intersect(rownames(X), rownames(Y))
      Y_proj <- Y
      
      #Calculate meaningful distance metric
      pr_X <- princomp(X)
      
      #We have a problme with d. When d is larger than the column of the data set we get an out of bound
      X_dist <- as.matrix(dist(pr_X$scores[,01:d], method = 'euclidean'))
      
      
      for (j in all_missing_proteins){
        
        protein_to_project <- j
        
        #Find the k closest neighbors that are also available in Y
        closest_proteins_index <- sort(X_dist[protein_to_project, all_common_proteins])[1:k]
        closest_proteins <- names(sort(X_dist[protein_to_project, all_common_proteins])[1:k])
        
        #Calculate the batch vector
        batch_vec <- X[closest_proteins,] - Y[closest_proteins,]
        
        #Project protein into Y batch space for each remaining neighbor
        projected_vec <- matrix(0, nrow = nrow(batch_vec), ncol = ncol(batch_vec))
        for (i in 1:nrow(batch_vec)){
          projected_vec[i,] <- X[protein_to_project,] - batch_vec[i,]
        }
        
        #Use the average projected vector
        average_proj <- matrix(colMeans(projected_vec), nrow = 1)
        
        # average_proj is a 1×#samples matrix
        rownames(average_proj) <- protein_to_project
        
        #Put the new protein into the batch corrected matrix
        Y_proj <- rbind(Y_proj, average_proj)
        
      }
      
      return(Y_proj)
      
    }
    
    # Iterate over the names of data_stes
    for (dataset_name in names(data_stes)) {
      
      X_count <- as.matrix(data_stes[[dataset_name]][["normalized"]])
      X_count[is.na(X_count)] <- 0
      
      # Extract the relevant pattern for the current dataset name
      pattern <- patterns_list[[dataset_name]]  # Match pattern with the dataset name
      
      # Use grep to extract columns based on the pattern
      Y_count <- combat_exprs[, grep(pattern, colnames(combat_exprs), ignore.case = TRUE, value = TRUE)]
      
      # Project the data using your custom function
      data_stes[[dataset_name]][["projection"]] <- as.data.frame(project_X_to_Y_proteins(X_count, Y_count))
      
    }
    
    ##Projection Matrix and Imputation
    all_genes <- lapply(data_stes, function(x){
      return(rownames(x$projection))
    })
    
    all_genes <- Reduce(union, all_genes)
    
    n <- names(data_stes)
    mat_names <- unlist(lapply(n, function(x) strsplit(x, "_")[[1]][1]))
    all_matrices <- lapply(data_stes, function(x){
      m <- x$projection[all_genes,]
      rownames(m) <- all_genes
      return(m)
    })
    
    for(i in 1:length(mat_names)){
      colnames(all_matrices[[i]]) <- paste0(mat_names[i],"-",colnames(all_matrices[[i]]))
    }
    
    all_matrices_proj <- Reduce(cbind, all_matrices)
    
    print("Projection Complete")
    
    #HeatMap
    # Replace NA with a specific value to identify them later
    all_matrices_mat_na <- all_matrices_proj %>% 
      mutate(across(everything(), ~replace(., is.na(.), -Inf)))
    
    # Gather the data into long format for ggplot
    long_data <- all_matrices_mat_na %>%
      rownames_to_column(var = "Row") %>%
      gather(key = "Column", value = "Value", -Row)
    
    # Create a new variable to identify the NA and zero values
    long_data <- long_data %>%
      mutate(
        ValueType = case_when(
          Value == -Inf ~ "NA",
          Value == 0 ~ "Zero",
          TRUE ~ "Other"
        )
      )
    
    # Define custom colors
    custom_colors <- c("NA" = "grey", "Zero" = "red", "Other" = "blue")
    
    # Create the heatmap
    heatmap_one <- ggplot(long_data, aes(x = Column, y = Row, fill = ValueType)) +
      geom_tile() +
      scale_fill_manual(values = custom_colors) +
      labs(title = "Heatmap of NA and Zero Values", fill = "Value Type") +
      theme_minimal() +
      theme(axis.text.x = element_blank())+
      theme(axis.text.y = element_blank())
    
    ggsave(
      filename = file.path(plots_dir, "Heatmap of NA.png"),
      plot     = heatmap_one,
      width    = 8,      # inches → 800 px at dpi=100
      height   = 6,      # inches → 600 px at dpi=100
      units    = "in",
      dpi      = 100
    )
    
    
    
    missing_percentage <- rowSums(is.na(all_matrices_proj)) / ncol(all_matrices_proj) * 100
    
    rows_to_impute <- missing_percentage <= 22
    matrix_to_impute <- as.matrix(all_matrices_proj[rows_to_impute, ])
    
    msnset_object <- MSnSet(exprs = matrix_to_impute)
    
    imputed_result <- imputeLCMD::impute.QRILC(exprs(msnset_object))
    
    imputed_matrix <- imputed_result[[1]]
    
    
    # Verify
    sum(is.na(imputed_matrix))   # Should be 0 if imputation was successful
    
    all_tissue_types <- c()
    for (l in names(data_stes)) {
      tissue_type <- data_stes[[l]][["meta"]][["tissue"]]
      all_tissue_types <- c(all_tissue_types, tissue_type)
    }
    all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))
    
    
    # Perform PCA - After Batch Correction
    pca_after <- prcomp(t(imputed_matrix), scale.=TRUE)
    pca_data_after <- data.frame(
      PC1 = pca_after$x[, 1], 
      PC2 = pca_after$x[, 2], 
      Batch = factor(batch_data$batch, levels = names(batch_names), labels = batch_names)
    )
    imputed_pca <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +  
      scale_color_manual(values = palette) +
      ggtitle("PCA Plot After Imputation 18%") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      ) +
      labs(color = "Batch")
    
    ggsave(
      filename = file.path(plots_dir, "PCA Plot After Imputation 18%.png"),
      plot     = imputed_pca,
      width    = 8,      # inches → 800 px at dpi=100
      height   = 6,      # inches → 600 px at dpi=100
      units    = "in",
      dpi      = 100
    )
    
    
    
    pca_data_after <- data.frame(
      PC1 = pca_after$x[, 1], 
      PC2 = pca_after$x[, 2], 
      Batch = all_tissue_types
    )
    tissue_impute <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +  
      scale_color_manual(values = palette) +  
      ggtitle("PCA Tissue 18% Imputation") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      ) +
      labs(color = "Tissues")
    ggsave(
      filename = file.path(plots_dir, "Tissue Plot After Imputation 18%.png"),
      plot     = tissue_impute,
      width    = 8,      # inches → 800 px at dpi=100
      height   = 6,      # inches → 600 px at dpi=100
      units    = "in",
      dpi      = 100
    )
    
    
    
    all_tissue_types <- c()
    for (l in names(data_stes)) {
      tissue_type <- data_stes[[l]][["meta"]][["diet"]]
      all_tissue_types <- c(all_tissue_types, tissue_type)
    }
    all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))
    
    pca_data_after <- data.frame(
      PC1 = pca_after$x[, 1], 
      PC2 = pca_after$x[, 2], 
      Batch = all_tissue_types
    )
    ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +  
      scale_color_manual(values = palette) +  
      ggtitle("PCA Diet 18% Imputation") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      ) +
      labs(color = "Diet")
    
    # Save final outputs (RDS)
    out_rds <- file.path(OUTPUT_DIR, paste0("imputed_matrix_", Sys.Date(), ".rds"))
    saveRDS(imputed_matrix, file = out_rds)
    message("Saved imputed matrix to: ", out_rds)
    
    # Save last plot PNG as well
    out_plot <- file.path(plot_dir, "last_plot.png")
    ggsave(filename = out_plot, plot = last_plot(), width = 10, height = 7, dpi = 100)
    message("Saved last_plot to: ", out_plot)
    
    # Return structured result for the app (matrix + plot_dir + small gallery)
    list(
      imputed_matrix = imputed_matrix,
      plot_dir = plot_dir,
      p_before = p_before,
      p_after = p_after,
      p_imputed = imputed_pca
    )
  }) # end processed_data eventReactive
  
  # This observer triggers when button is pressed (reads the cached eventReactive result)
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
    # store results
    results(pd)
    output$analysis_status <- renderText("✅ Analysis finished.")
    
    # collect PNGs from plot_dir (if any) and publish to UI
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
  
  # Plots gallery
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
      paste0("processed_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      dat <- results(); req(!is.null(dat))
      if (is.data.frame(dat)) {
        write.csv(dat, file, row.names = FALSE)
      } else if (is.list(dat) && !is.null(dat$imputed_matrix)) {
        write.csv(as.data.frame(dat$imputed_matrix), file, row.names = FALSE)
      } else {
        writeLines("Output is not a data frame.", file)
      }
    }
  )
} # end server

shinyApp(ui = ui, server = server)

