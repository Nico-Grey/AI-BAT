library(shiny)
library(DT)      # For displaying data tables
library(ggplot2) # For plotting
library(png)     # For reading PNG images
library(grid)    # For displaying plots
library(shinyjs)

ui <- navbarPage("AI-BATS",
                 shinyjs::useShinyjs(),
                 tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"))
                 ,
                 # 🏠 HOME TAB
                 tabPanel("Home",
                          fluidRow(
                            column(2, tags$img(src = "logo.JFIF", width = "100px")),
                            column(10, h3("Welcome to AI-BATS"))
                          ),
                          p("AI-Bats Placeholder text")
                 ),
                 
                 # 📂 UPLOAD DATA TAB
                 tabPanel("Upload Data",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("file", "Upload Data File", accept = c(".csv", ".RDS", ".RData")),
                              actionButton("run", "Run Processing")
                            ),
                            mainPanel(
                              h4("Instructions"),
                              p("Upload your dataset in .RDS or .RData format and click 'Run Processing'.")
                            )
                          )
                 ),
                 
                 # 📊 RESULTS TAB (PROCESSED DATA + IMAGES)
                 tabPanel("Results",
                          fluidRow(
                            column(6, 
                                   h4("Processed Data"),
                                   DTOutput("data_table")
                            ),
                            column(6,  
                                   h4("Before Processing"),
                                   tags$img(src = "before.png", width = "80%", style = "max-width: 250px; height: auto;"),
                                   h4("After Processing"),
                                   tags$img(src = "after.png", width = "80%", style = "max-width: 250px; height: auto;")
                            )
                          ),
                          plotOutput("plot")  # Plot at the bottom
                 ),
                 
                 # ⬇️ DOWNLOAD TAB
                 tabPanel("Download",
                          sidebarLayout(
                            sidebarPanel(
                              downloadButton("download", "Download Processed Data")
                            ),
                            mainPanel(
                              h4("Download Your Results"),
                              p("Click the button to download the processed dataset in .RDS format.")
                            )
                          )
                 ),
                 
                 # ⚙️ SETTINGS TAB
                 tabPanel("Settings",
                          sidebarLayout(
                            sidebarPanel(
                              checkboxInput("dark_mode", "Enable Dark Mode", FALSE),
                              sliderInput("image_size", "Image Size", min = 50, max = 400, value = 250)
                            ),
                            mainPanel(
                              h4("Customize Your Experience"),
                              p("Toggle settings like dark mode and image size.")
                            )
                          )
                 ),
                 
                 # ℹ️ ABOUT TAB
                 tabPanel("About",
                          h3("About AI-BATS"),
                          p("AI-BATS is an advanced data processing tool designed to automate data cleaning and visualization."),
                          p("")
                 )
)


server <- function(input, output, session) {
  ### DARK MODE FUNCTIONALITY
  observe({
    if (input$dark_mode) {
      shinyjs::addClass(selector = "body", class = "dark-mode")  # Apply dark mode
    } else {
      shinyjs::removeClass(selector = "body", class = "dark-mode")  # Remove dark mode
    }
  })
  
  results <- reactiveVal(NULL)  # Store processed data
  output_plot <- reactiveVal(NULL)  # Store processed plot file path
  
  # Reactive function to handle different file types
  uploaded_data <- reactive({
    req(input$file)  # Ensure a file is uploaded
    
    ext <- tools::file_ext(input$file$name)  # Get file extension
    
    # Read file based on its type
    data <- switch(ext,
                   "csv" = read.csv(input$file$datapath, stringsAsFactors = FALSE),
                   "RDS" = readRDS(input$file$datapath),
                   "RData" = {
                     temp_env <- new.env()
                     load(input$file$datapath, envir = temp_env)
                     get(ls(temp_env)[1], envir = temp_env)  # Return first object found
                   },
                   stop("Invalid file type")
    )
    
    validate(need(!is.null(data), "Failed to load data. Check file format."))  
    return(data)
  })
  
  observeEvent(input$run, {
    req(input$file)  # Ensure a file is uploaded
    
    # Convert data to RDS for Docker processing
    input_data <- uploaded_data()  # Get uploaded data
    input_path <- tempfile(fileext = ".RDS")  # Temporary input file
    saveRDS(input_data, input_path)  # Save as RDS
    
    # Define temporary output paths
    output_data_path <- tempfile(fileext = ".RDS")  # Temporary output file
    output_plot_path <- tempfile(fileext = ".png")  # Temporary plot file
    
    # Run script inside Docker
    docker_command <- sprintf(
      "docker run --rm -v %s:/app/input.RDS -v %s:/app/output.RDS -v %s:/app/output.png my_r_shiny_app Rscript /app/preprocess.R /app/input.RDS /app/output.RDS /app/output.png",
      input_path, output_data_path, output_plot_path
    )
    
    system(docker_command)  # Run the command
    
    # Load results back into R
    processed_data <- readRDS(output_data_path)
    results(processed_data)  # Store in reactive value
    output_plot(output_plot_path)  # Store plot file path
  })
  
  # Show the processed data as a table
  output$data_table <- renderDT({
    req(results())
    datatable(results())
  })
  
  # Show the generated plot
  output$plot <- renderPlot({
    req(output_plot())
    img <- readPNG(output_plot())
    grid::grid.raster(img)
  })
  
  # Provide a download option
  output$download <- downloadHandler(
    filename = function() { "processed_data.RDS" },
    content = function(file) {
      file.copy(results(), file)
    }
  )
}

  
 



shinyApp(ui, server)
