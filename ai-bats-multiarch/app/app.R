library(shiny)

# Define UI for the application
ui <- fluidPage(
    titlePanel("AI-BATS: Artificial Intelligence-driven Brown Adipose Tissue Scoring"),
    
    sidebarLayout(
        sidebarPanel(
            # Input elements can be added here
            h3("Input Parameters"),
            # Example input: file upload
            fileInput("datafile", "Upload Dataset", accept = ".rds"),
            actionButton("run_analysis", "Run Analysis")
        ),
        
        mainPanel(
            # Output elements can be added here
            h3("Results"),
            verbatimTextOutput("results")
        )
    )
)

# Define server logic
server <- function(input, output) {
    observeEvent(input$run_analysis, {
        req(input$datafile)
        
        # Load the dataset
        dataset <- readRDS(input$datafile$datapath)
        
        # Call the main analysis pipeline
        source("BATpipeline.R")
        results <- run_analysis(dataset)  # Assuming run_analysis is a function in BATpipeline.R
        
        # Output the results
        output$results <- renderPrint({
            results
        })
    })
}

# Run the application 
shinyApp(ui = ui, server = server)