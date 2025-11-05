#!/bin/bash
#
# Build and push Docker containers for LongReads2 pipeline
#
# Usage:
#   ./build_containers.sh [--push] [--registry YOUR_REGISTRY]
#
# Options:
#   --push      Push images to Docker registry after building
#   --registry  Specify custom registry (default: longreads2)

set -euo pipefail

# Default values
PUSH=false
REGISTRY="longreads2"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --push)
            PUSH=true
            shift
            ;;
        --registry)
            REGISTRY="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Building LongReads2 Docker containers${NC}"
echo "Registry: $REGISTRY"
echo ""

# Build containers
containers=("ncbi-tools" "sequence-analysis" "visualization")

for container in "${containers[@]}"; do
    echo -e "${GREEN}Building ${container}...${NC}"
    docker build -t ${REGISTRY}/${container}:latest ${container}/
    echo ""
done

echo -e "${GREEN}All containers built successfully!${NC}"
echo ""

# List built images
echo "Built images:"
docker images | grep "$REGISTRY"
echo ""

# Push if requested
if [ "$PUSH" = true ]; then
    echo -e "${BLUE}Pushing containers to registry...${NC}"
    for container in "${containers[@]}"; do
        echo -e "${GREEN}Pushing ${container}...${NC}"
        docker push ${REGISTRY}/${container}:latest
    done
    echo -e "${GREEN}All containers pushed successfully!${NC}"
fi

echo ""
echo "To run a test:"
echo "  docker run --rm -it ${REGISTRY}/ncbi-tools:latest datasets --version"
echo "  docker run --rm -it ${REGISTRY}/sequence-analysis:latest barrnap --version"
echo "  docker run --rm -it ${REGISTRY}/visualization:latest R --version"
