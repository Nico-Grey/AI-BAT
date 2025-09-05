#!/bin/bash
set -euo pipefail

IMAGE_NAME="lst93/aibats:app-test"


# Create or reuse builder
if ! docker buildx inspect multiarch-builder >/dev/null 2>&1; then
  docker buildx create --name multiarch-builder --driver docker-container --use
  docker buildx inspect --bootstrap
else
  docker buildx use multiarch-builder
fi

# (Usually not needed on Docker Desktop; keep if required)
# docker run --rm --privileged multiarch/qemu-user-static --reset -p yes

# Prepare local cache directories
mkdir -p .buildx-cache .buildx-cache-new

docker buildx build \
  --platform linux/amd64 \
  -t "${IMAGE_NAME}" \
  --push \
  --build-arg BUILDKIT_INLINE_CACHE=1 \
  --cache-from type=local,src=.buildx-cache \
  --cache-to type=local,dest=.buildx-cache-new,mode=max \
  . 

# Atomically replace old cache with new
rm -rf .buildx-cache
mv .buildx-cache-new .buildx-cache

echo "Done. Next build will reuse cache. To inspect:"
echo "  docker buildx imagetools inspect ${IMAGE_NAME}"