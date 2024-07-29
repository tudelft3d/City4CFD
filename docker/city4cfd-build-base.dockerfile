FROM ubuntu:noble
LABEL org.opencontainers.image.authors="Ivan Paden <i.paden@tudelft.nl>"
LABEL org.opencontainers.image.source="https://github.com/tudelft3d/city4cfd"
LABEL org.opencontainers.image.vendor="Tudelft3D"
LABEL org.opencontainers.image.title="City4CFD build base"
LABEL org.opencontainers.image.description="Base image for building City4CFD"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.url="https://github.com/tudelft3d/city4cfd"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential    \
    cmake              \
    libboost-all-dev   \
    libcgal-dev        \
    libeigen3-dev      \
    libomp-dev         \
    libgdal-dev
