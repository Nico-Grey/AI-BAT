library(shiny)
library(DT)      # For displaying data tables
library(ggplot2) # For plotting
library(png)     # For reading PNG images
library(grid)    # For displaying plots

ui <- fluidPage(
  # Add a title with a logo
  fluidRow(
    column(1, 
           tags$img(src = "logo.JFIF", height = "250px")  # Adjust path & size
    ),
    column(11, 
           titlePanel("AI-BATS")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload R Data File", accept = c(".RDS", ".RData")),
      actionButton("run", "Run Processing"),
      downloadButton("download", "Download Processed Data")
    ),
    
    mainPanel(
      fluidRow(
        column(6,  # Right: Before/After Images (closer to the data table)
               h4("Before Processing"),
               tags$img(src = "before.png", width = "100%", style = "max-width: 250px; height: auto; display: block; margin-left: auto; margin-right: auto;")
        ),
        column(6, h4("After Processing"),
               tags$img(src = "after.png", width = "100%", style = "max-width: 250px; height: auto; display: block; margin-left: auto; margin-right: auto;")
        )
      ),
      plotOutput("plot")  # Keep the plot below if needed
    )
    
    
  )
)

server <- function(input, output, session) {
  
  results <- reactiveVal(NULL)  # Store processed data
  output_plot <- reactiveVal(NULL) # Store processed plot file path
  
  observeEvent(input$run, {
    req(input$file)  # Ensure a file is uploaded
    
    # Define temporary paths for mounting files
    input_path <- input$file$datapath
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
