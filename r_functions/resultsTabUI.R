resultsTabUI <- function() {
  tabsetPanel(
    id = "results_tabs",
    tabPanel(
      title = "Processed Data",
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
            width = 4,
            h3("Plots"),
            p("Interactive gallery of output plots."),
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("a) Plot 1", style = "text-align:center; font-size:14px;"),
              imageOutput("plot1", width = "300px", height = "250px")
            )
          ),
          column(
            width = 4,
            h3("Plots"),
            p("Interactive gallery of output plots."),
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("b) Plot 2", style = "text-align:center; font-size:14px;"),
              imageOutput("plot2", width = "300px", height = "250px")
            )
          ),
          column(
            width = 4,
            h3("Plots"),
            p("Interactive gallery of output plots."),
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("c) Plot 3", style = "text-align:center; font-size:14px;"),
              imageOutput("plot3", width = "300px", height = "250px")
            )
          )
        )
      )
    ),
    tabPanel(
      title = "Final Results",
      fluidPage(
        fluidRow(
          column(
            width = 12,
            h3("Final Results Summary"),
            p("This section displays the final processed outcomes, analysis metrics, or combined summaries."),
            verbatimTextOutput("final_summary"),
            br(),
            h4("Final Results Table"),
            DTOutput("final_results_table")
          )
        ),
        fluidRow(
          column(
            width = 12,
            h4("Final Figures"),
            div(
              class = "plot-gallery",
              # placeholder_image_card("Final Plot 1", "placeholder_plot1.png"),
              # placeholder_image_card("Final Plot 2", "placeholder_plot2.png")
            )
          )
        )
      )
    )
  )
}