# AI-BATS Multi-Architecture Project

This project provides a containerized environment for the AI-BATS (Artificial Intelligence-driven Brown Adipose Tissue Scoring) application, which includes R, Shiny, Python, and all necessary dependencies for protein expression analysis and machine learning workflows.

## Project Structure

- **app/**: Contains the main application files.
  - `app.R`: Shiny web application interface.
  - `BATpipeline.R`: Main analysis pipeline script.
  - `general_functions.R`: Utility functions for analysis.
  - **data/**: Directory containing datasets.
    - `dataset.rds`: Example/test dataset for analysis.

- **docker/**: Contains Docker-related files.
  - `Dockerfile`: Configuration for building the Docker image.
  - **build/**: Scripts for building multi-architecture images.
    - `buildx.sh`: Script to build Docker images for x86 and ARM architectures.
    - `qemu-setup.sh`: Script to set up QEMU for emulating different architectures.
  - `README.md`: Documentation for Docker setup.

- `environment.yml`: Specifies the Python environment and dependencies.

## Building Multi-Architecture Images

To build Docker images for both x86 and ARM architectures, follow these steps:

1. **Set Up QEMU**: Run the `qemu-setup.sh` script to configure QEMU for emulation.
   ```bash
   cd docker/build
   ./qemu-setup.sh
   ```

2. **Build the Images**: Use the `buildx.sh` script to build the images.
   ```bash
   ./buildx.sh
   ```

## Running the Application

After building the Docker image, you can run the application using the following command:
```bash
docker run -p 3838:3838 my_r_shiny_app
```

## Usage

For detailed usage instructions and additional information, refer to the documentation in the respective directories.