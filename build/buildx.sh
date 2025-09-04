#!/bin/bash

# Set the Docker image name
IMAGE_NAME="lst93/aibats:app"

# Create a builder instance for multi-architecture builds
docker buildx create --name multiarch-builder --use

# Set up QEMU for emulating different architectures
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes

# Build the Docker image for multiple architectures
docker buildx build --platform linux/amd64,linux/arm64 -t $IMAGE_NAME-latest --push .

# Clean up the builder instance
docker buildx rm multiarch-builder