library(shiny)
library(bslib)
library(DT)

# Define theme with light and dark blue color scheme
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

# UI for Data Input tab
inputTabUI <- function() {
  fluidPage(
    fluidRow(
      column(
        width = 4,
        h3("Upload Data"),
        p("Select your input data file below (CSV or RDS format)."),
        fileInput("data_file", "Choose File", 
                  accept = c(".csv", ".rds"),
                  buttonLabel = "Browse...", placeholder = "No file selected"),
        selectInput("plot_color", "Choose Point Color:",
                    choices = c("blue", "red", "green", "purple", "black"),
                    selected = "blue"),
        sliderInput("plot_size", "Point Size:",
                    min = 1, max = 10, value = 3),
        actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary")
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

# UI for Results tab
resultsTabUI <- function() {
  fluidPage(
    fluidRow(
      column(
        width = 12,
        h3("Processed Data"),
        p("Below is the table of processed data. Download it using the button below."),
        div(class = "data-table-container",
            DTOutput("results_table")
        ),
        downloadButton("download_data", "Download Data", class = "btn-secondary")
      )
    ),
    fluidRow(
      column(
        width = 12,
        h3("Plots"),
        p("Interactive gallery of output plots. Resize the window for responsive layout."),
        div(class = "plot-gallery",
            uiOutput("plots_gallery")
        )
      )
    )
  )
}

# UI for Examples & Tutorial tab
examplesTabUI <- function() {
  fluidPage(
    h2("Examples & Tutorial"),
    p("This section provides example usage, sample data, and instructions."),
    h3("Step-by-Step Instructions"),
    tags$ul(
      tags$li("Upload a CSV or RDS file using the Data Input tab."),
      tags$li("Click 'Run Analysis' to process the data."),
      tags$li("Review the results in the Results tab."),
      tags$li("Download processed data as needed.")
    ),
    h3("Example Input Data"),
    p("Replace this with a description or snippet of your example input data format."),
    h3("Sample Output"),
    fluidRow(
      column(
        width = 6,
        h4("Sample Plot 1"),
        p("Placeholder for an example output plot or figure.")
      ),
      column(
        width = 6,
        h4("Sample Plot 2"),
        p("Placeholder for an additional example output plot or figure.")
      )
    )
  )
}

# Main UI
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

# placeholders for results & plots


# Server logic
server <- function(input, output, session) {
  
  results     <- reactiveVal(NULL)
  output_plots <- reactiveVal(character(0))
  
  options(shiny.maxRequestSize = 150 * 1024^2)

  # Reactive data input
  data_input <- reactive({
    req(input$data_file)
    ext <- tolower(tools::file_ext(input$data_file$name))
    switch(ext,
           csv = read.csv(input$data_file$datapath, header = TRUE, stringsAsFactors = FALSE),
           rds = readRDS(input$data_file$datapath),
           stop("Unsupported file type.")
    )
  })

  # Preview table
  output$preview_table <- renderDT({
    dat <- data_input(); req(is.data.frame(dat))
    datatable(head(dat, 5), options = list(dom = 't', paging = FALSE))
  })
  
  # Preview non-tabular objects
  output$non_table_preview <- renderPrint({
    obj <- data_input()
    if (!is.data.frame(obj)) str(obj)
  })
  
  # Theme toggle
  observeEvent(input$theme_choice, {
    newt <- if (input$theme_choice == "Light") {
      bs_theme_update(my_theme, bootswatch = "cerulean", primary = "#005f99")
    } else {
      bs_theme_update(my_theme, bootswatch = "darkly",   primary = "#005f99")
    }
    (bs_theme(version = 4))
  })
  
  # Reactive holders for results and plots
  # Reactive holders for results & plots
  results     <- reactiveVal(NULL)
  output_plots <- reactiveVal(character(0))
  
  #Test for directory
  message("R’s tempdir() is: ", tempdir())
  message("Contents of tempdir():")
  print(list.files(tempdir(), full.names = TRUE))
  
  
  # Main “Run analysis” handler
  observeEvent(input$run_analysis, {
    req(data_input())
    
    # 1) Snapshot inputs & params
    tmp_in     <- tempfile(fileext = ".RDS")
    tmp_params <- tempfile(fileext = ".RDS")
    saveRDS(data_input(), tmp_in)
    saveRDS(list(color = input$plot_color, size = input$plot_size), tmp_params)
    
    # 2) Prepare output locations
    tmp_out   <- tempfile(fileext = ".RDS")
    plots_dir <- tempfile(); dir.create(plots_dir)
    
    # 3) Normalize paths & build Docker cmd
    #in_path  <- gsub("\\\\", "/", tmp_in)
    #out_path <- gsub("\\\\", "/", tmp_out)
    #plots_m  <- gsub("\\\\", "/", plots_dir)
    #param_p  <- gsub("\\\\", "/", tmp_params)
    
    
    docker_cmd <- paste0(
      'docker run --rm ',
      '-v "', tmp_in,     '":/app/input.RDS ',
      '-v "', tmp_out,    '":/app/output.RDS ',
      '-v "', plots_dir,  '":/app/plots ',
      '-v "', tmp_params, '":/app/params.RDS ',
      'my_r_shiny_app ',
      'Rscript /app/BATpipeline.R /app/input.RDS /app/output.RDS /app/plots /app/params.RDS'
    )
    
    
    # 4) Execute and load results
    system(docker_cmd)
    results(readRDS(tmp_out))
    output_plots(list.files(plots_dir, pattern = "\\.png$", full.names = TRUE))
  }, ignoreInit = TRUE)
  
  # Show processed data in Results tab
  output$results_table <- renderDT({
    dat <- results(); req(is.data.frame(dat))
    datatable(dat, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  # Serve plots via renderImage inside renderUI
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
  
  # Download handler for processed data
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("processed_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      dat <- results(); req(!is.null(dat))
      if (is.data.frame(dat)) {
        write.csv(dat, file, row.names = FALSE)
      } else {
        writeLines("Output is not a data frame.", file)
      }
    }
  )
  
}  # <-- close server()

# Launch the app
shinyApp(ui = ui, server = server)
