# AI-BATS Shiny Application

## Overview

AI-BATS (Artificial Intelligence-driven Brown Adipose Tissue Scoring) is an interactive Shiny web application designed for analyzing protein expression data related to brown and white adipose tissue. The application provides a user-friendly interface for uploading data, running computational analysis pipelines, and visualizing results.

## Features

### 🎯 Core Functionality
- **Data Upload**: Support for CSV and RDS file formats
- **Interactive Analysis**: One-click analysis execution with customizable parameters
- **Real-time Visualization**: Dynamic plot generation and display
- **Data Export**: Download processed results in CSV format
- **Responsive Design**: Modern, mobile-friendly interface

### 📊 Analysis Capabilities
- Protein expression data processing
- Brown/white adipose tissue classification
- Machine learning-based scoring algorithms
- Dimensionality reduction and visualization
- Batch correction and normalization

### 🎨 User Interface
- **Dark/Light Theme Toggle**: Switch between dark and light modes
- **Tabbed Navigation**: Organized workflow across multiple tabs
- **Interactive Gallery**: Responsive plot display with automatic layout
- **Data Preview**: Real-time preview of uploaded data
- **Progress Feedback**: Visual indicators for long-running processes

## Installation

### Prerequisites
- R (version 4.0 or higher)
- Docker (for containerized analysis pipeline)
- Required R packages (automatically installed):
  - `shiny`
  - `bslib`
  - `DT`

### Setup
1. Clone the repository
2. Install required R packages:
```r
install.packages(c("shiny", "bslib", "DT"))
```
3. Build the Docker container for the analysis pipeline:
```bash
docker build -t my_r_shiny_app .
```

## Usage

### Starting the Application
```r
# From R console
shiny::runApp("app.R")
```

### Workflow

#### 1. Data Input Tab
- **Upload Data**: Select your protein expression file (CSV or RDS format)
- **Customize Visualization**: 
  - Choose point colors for plots
  - Adjust point sizes (1-10 scale)
- **Run Analysis**: Click the "Run Analysis" button to process your data
- **Data Preview**: View the first 5 rows of your uploaded data

#### 2. Results Tab
- **Processed Data Table**: Interactive table with processed results
- **Plot Gallery**: Automatically generated visualizations
- **Download**: Export processed data as CSV

#### 3. Examples & Tutorial Tab
- Step-by-step instructions
- Example data formats
- Sample outputs and interpretations

#### 4. Settings
- **Appearance**: Toggle between dark and light themes

### Data Format Requirements

#### Input Data Structure
Your input file should contain:
- **Rows**: Samples/observations
- **Columns**: Protein identifiers or gene symbols
- **Values**: Expression levels or intensity measurements

#### Supported File Formats
- **CSV**: Comma-separated values with headers
- **RDS**: R data serialization format

#### Example Data Structure
```
Sample_ID | Protein_A | Protein_B | Protein_C | ...
----------|-----------|-----------|-----------|----
Sample_1  |    245.3  |    18.7   |    892.1  | ...
Sample_2  |    312.8  |    24.1   |    756.4  | ...
...       |    ...    |    ...    |    ...    | ...
```

## Technical Architecture

### Frontend (UI)
- **Framework**: Shiny with Bootstrap 4 styling
- **Theme System**: Custom dark/light theme implementation
- **Layout**: Responsive grid system with CSS Grid for plot galleries
- **Components**: Interactive data tables, file inputs, plot outputs

### Backend (Server)
- **Data Processing**: Reactive programming model
- **File Handling**: Support for multiple data formats
- **Analysis Pipeline**: Docker-containerized R processing
- **Memory Management**: Optimized for large datasets (up to 150MB)

### Integration
- **Docker Integration**: Seamless container execution for analysis
- **File System**: Temporary file management for secure processing
- **Path Handling**: Cross-platform path normalization

## Configuration

### Memory Limits
```r
options(shiny.maxRequestSize = 150 * 1024^2)  # 150MB limit
```

### Docker Mount Points
The application automatically mounts:
- Input data: `/app/input.RDS`
- Output results: `/app/output.RDS`
- Plot directory: `/app/plots`
- Parameters: `/app/params.RDS`

## Troubleshooting

### Common Issues

#### File Upload Problems
- **Issue**: File not uploading
- **Solution**: Check file size (must be < 150MB) and format (CSV/RDS only)

#### Analysis Not Running
- **Issue**: "Run Analysis" button not working
- **Solution**: Ensure Docker is running and container is built

#### Plots Not Displaying
- **Issue**: Empty plot gallery
- **Solution**: Check that analysis completed successfully and plot files were generated

#### Theme Not Switching
- **Issue**: Dark/light theme toggle not working
- **Solution**: Refresh the application and try again

### Error Messages
- **"Unsupported file type"**: Use only CSV or RDS formats
- **"Docker command failed"**: Check Docker installation and container build
- **"Output is not a data frame"**: Analysis pipeline may have failed

## Performance Optimization

### Large Datasets
- Use RDS format for faster loading
- Consider data preprocessing before upload
- Monitor memory usage during analysis

### Plot Performance
- Plots are cached for faster display
- Use appropriate point sizes for large datasets
- Gallery layout automatically adjusts for screen size

## Development

### Code Structure
```
app.R
├── UI Definition
│   ├── Theme Configuration
│   ├── Tab Layouts
│   └── CSS Styling
└── Server Logic
    ├── Data Input Handling
    ├── Analysis Pipeline
    ├── Plot Generation
    └── Download Handlers
```

### Extending Functionality
To add new features:
1. Modify the UI layout in the appropriate tab function
2. Add server logic for new reactive elements
3. Update Docker pipeline if needed
4. Test with sample data

## Support

For issues, questions, or contributions:
1. Check the troubleshooting section above
2. Review example data in the Examples & Tutorial tab
3. Ensure all prerequisites are properly installed
4. Verify Docker container is functioning correctly

## Version Information

- **Current Version**: 1.0
- **R Version**: 4.0+
- **Shiny Version**: 1.7+
- **Bootstrap Version**: 4.x
