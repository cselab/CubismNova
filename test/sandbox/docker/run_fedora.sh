set -e
docker run -it --mount type=bind,source="$(pwd -P)",target=/test --rm cselab/fedora_mpich:latest /bin/bash
