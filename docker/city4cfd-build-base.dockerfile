FROM ubuntu:noble
LABEL org.opencontainers.image.authors="Ivan Paden <i.paden@tudelft.nl>"
LABEL org.opencontainers.image.source="https://github.com/tudelft3d/city4cfd"
LABEL org.opencontainers.image.vendor="Tudelft3D"
LABEL org.opencontainers.image.title="City4CFD build base"
LABEL org.opencontainers.image.description="Base image for building City4CFD"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.url="https://github.com/tudelft3d/city4cfd"

# Install required packages and clean up
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential    \
    cmake              \
    libboost-all-dev   \
    libgmp-dev         \
    libmpfr-dev        \
    libeigen3-dev      \
    libomp-dev         \
    libgdal-dev        \
    wget               \
    ca-certificates    \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*

# Download and extract CGAL
RUN wget https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1-library.tar.xz \
    && tar -xf CGAL-6.0.1-library.tar.xz -C /opt \
    && rm CGAL-6.0.1-library.tar.xz

# Set environment variable for CGAL
ENV CGAL_DIR=/opt/CGAL-6.0.1
