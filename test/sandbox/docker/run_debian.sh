set -e
docker run -it --mount type=bind,source="$(pwd -P)",target=/test --rm cselab/debian_buster_mpich:latest /bin/bash
