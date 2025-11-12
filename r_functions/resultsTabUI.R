resultsTabUI <- function() {
  tabsetPanel(
    id = "results_tabs",
    
    # ---- PROCESSED DATA TAB ----
    tabPanel(
      title = "Processed Data",
      fluidPage(
        # Section header
        fluidRow(
          column(
            width = 12,
            h3("Processed Data"),
            downloadButton("download_data", "Download Data", class = "btn-secondary"),
            br(), br()
          )
        ),
        
        # ---- FIRST ROW OF PLOTS ----
        fluidRow(
          column(
            width = 4,
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("Normalization per dataset", style = "text-align:center; font-size:14px;"),
              imageOutput("plot1", width = "100%", height = "250px")
            )
          ),
          column(
            width = 4,
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("PCA after normalization", style = "text-align:center; font-size:14px;"),
              imageOutput("plot2", width = "100%", height = "250px")
            )
          ),
          column(
            width = 4,
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("PCA after COMBAT", style = "text-align:center; font-size:14px;"),
              imageOutput("plot3", width = "100%", height = "250px")
            )
          )
        ),
        
        # ---- SECOND ROW OF PLOTS ----
        fluidRow(
          column(
            width = 4,
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("PCA after projection", style = "text-align:center; font-size:14px;"),
              imageOutput("plot4", width = "100%", height = "250px")
            )
          ),
          column(
            width = 4,
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("PCA after imputation", style = "text-align:center; font-size:14px;"),
              imageOutput("plot5", width = "100%", height = "250px")
            )
          ),
          column(
            width = 4,
            div(
              class = "card shadow-sm mb-3 p-2",
              h4("PCA after imputation by tissue", style = "text-align:center; font-size:14px;"),
              imageOutput("plot6", width = "100%", height = "250px")
            )
          )
        )
      )
    ),
    
    # ---- FINAL RESULTS TAB ----
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
              class = "plot-gallery"
              # placeholder_image_card("Final Plot 1", "placeholder_plot1.png"),
              # placeholder_image_card("Final Plot 2", "placeholder_plot2.png")
            )
          )
        )
      )
    )
  )
}
