#!/bin/bash

# This script sets up QEMU for multi-architecture builds using Docker Buildx.
# It enables the emulation of different CPU architectures, allowing the building
# of images for both x86 and ARM architectures.

set -e

# Install QEMU if not already installed
if ! command -v qemu-arm-static &> /dev/null; then
    echo "Installing QEMU..."
    apt-get update && apt-get install -y qemu-user-static
fi

# Register QEMU for the supported architectures
echo "Registering QEMU for multi-architecture support..."
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes

echo "QEMU setup complete. You can now build multi-architecture images using Docker Buildx."