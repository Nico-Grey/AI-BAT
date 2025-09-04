# AI-BATS Multi-Architecture Docker Setup

This directory contains the Docker configuration and scripts for building and running the AI-BATS (Artificial Intelligence-driven Brown Adipose Tissue Scoring) application on multiple architectures, including x86 and ARM.

## Overview

The AI-BATS application is containerized using Docker, allowing for easy deployment and scalability. This setup leverages Docker Buildx and QEMU to facilitate the building of images that can run on different hardware architectures.

## Directory Structure

- **Dockerfile**: Defines the Docker image configuration for the AI-BATS application, including the installation of R, Shiny, Python, and all necessary dependencies.
- **build/**: Contains scripts for building multi-architecture images.
  - **buildx.sh**: Script to build Docker images for multiple architectures using Docker Buildx.
  - **qemu-setup.sh**: Script to set up QEMU for emulating different architectures, enabling the building of multi-architecture images.

## Building Multi-Architecture Images

To build the Docker images for both x86 and ARM architectures, follow these steps:

1. **Install Docker**: Ensure you have Docker installed on your machine. You may also need to enable experimental features to use Buildx.

2. **Set Up QEMU**: Run the `qemu-setup.sh` script to set up QEMU for cross-platform builds. This will allow Docker to emulate different architectures.

   ```bash
   ./build/qemu-setup.sh
   ```

3. **Build the Images**: Use the `buildx.sh` script to build the Docker images for the desired architectures.

   ```bash
   ./build/buildx.sh
   ```

4. **Run the Application**: Once the images are built, you can run the AI-BATS application using Docker.

   ```bash
   docker run -p 3838:3838 <image_name>
   ```

Replace `<image_name>` with the name of the built image.

## Additional Information

For more details on the AI-BATS application, refer to the main `README.md` file located in the root of the project directory.