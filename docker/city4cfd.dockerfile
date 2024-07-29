FROM tudelft3d/city4cfd:build-base AS builder

ARG JOBS

#
# Install City4CFD
#
COPY . /tmp

RUN cd /tmp && \
    mkdir build && \
    cd build && \
    cmake \
#        -DCGAL_DIR=/usr/local \
#        -DCMAKE_BUILD_TYPE=Release \
#        -DBoost_NO_BOOST_CMAKE=TRUE \
#        -DBoost_NO_SYSTEM_PATHS=TRUE \
#        -DBOOST_ROOT=/usr/local \
        .. && \
    make -j $JOBS && \
    make install && \
    cd ~ && \
    rm -rf /tmp/* && \
    rm -rf /user/local/man

RUN city4cfd --version

# removing unnecessary headers
RUN rm -rf /usr/local/include

RUN mkdir /data && \
    chown 1001 /data && \
    chgrp 0 /data && \
    chmod g=u /data && \
    chgrp 0 /etc/passwd && \
    chmod g=u /etc/passwd

#
# Export the dependencies
#
RUN mkdir /export
COPY docker/strip-docker-image-export /tmp
RUN bash /tmp/strip-docker-image-export \
    -v \
    -d /export \
    -f /bin/bash \
    -f /usr/bin/awk \
    -f /usr/bin/id \
    -f /etc/passwd \
    -f /bin/ls \
    -f /data \
    -f /usr/share/proj/proj.db \
    -f /usr/local/bin/city4cfd_pcprep \
    -f /usr/local/bin/city4cfd_las2las \
    -f /usr/local/bin/city4cfd

#
# Create City4CFD image
#
FROM scratch AS exe
ARG VERSION
LABEL org.opencontainers.image.authors="Ivan Paden <i.paden@tudelft.nl>"
LABEL org.opencontainers.image.source="https://github.com/tudelft3d/city4cfd"
LABEL org.opencontainers.image.vendor="Tudelft3D"
LABEL org.opencontainers.image.title="City4CFD"
LABEL org.opencontainers.image.description="City4CFD image"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.url="https://github.com/tudelft3d/city4cfd"
LABEL org.opencontainers.image.version=$VERSION

COPY --from=builder /export/ /

WORKDIR /data

#ENTRYPOINT ["city4cfd"]
#CMD ["--help"]
