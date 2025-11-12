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
      # column(width = 6, placeholder_image_card("Sample Plot 1", "placeholder_plot1.png")),
      # column(width = 6, placeholder_image_card("Sample Plot 2", "placeholder_plot2.png"))
    )
  )
}