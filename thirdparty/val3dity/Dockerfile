FROM alpine:3.10
LABEL maintainer.email="b.dukai@tudelft.nl" maintainer.name="Balázs Dukai"

ENV CXX="g++ -std=c++98"

# Install geos
RUN apk --update add --virtual .geos-deps \
        which \
        make \
        gcc \
        g++ \
        file \
        git \
        autoconf \
        automake \
        libtool && \
    cd /tmp && \
    git clone https://git.osgeo.org/gitea/geos/geos.git geos && \
    cd geos && \
    git checkout 3.7.1 && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd / && \
    apk del .geos-deps && \
    rm -rf /tmp/* && \
    rm -rf /user/local/man

# install CGAL
RUN apk --update add --virtual .cgal-deps \
        which \
        make \
        cmake \
        gcc \
        g++ \
        eigen-dev \
        boost-dev \
        gmp-dev \
        mpfr-dev \
        zlib-dev && \
    cd /tmp && \
    wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.2/CGAL-5.0.2.tar.xz && \
    tar xf CGAL-5.0.2.tar.xz && \
    rm -f CGAL-5.0.2.tar.xz && \
    cd CGAL-5.0.2 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DWITH_Eigen3=ON .. && \
    make && \
    make install && \
    cd / && \
    apk del .cgal-deps && \
    rm -rf /tmp/* && \
    rm -rf /user/local/man

# install val3dity
COPY . /tmp
RUN apk --update add --virtual .val3dity-deps \
        make \
        cmake \
        gcc \
        g++ \
        boost-dev \
        eigen-dev \
        gmp-dev \
        mpfr-dev \
        python3 \
        pytest \
        py-yaml \
        py-lxml \
        py3-setuptools && \
    cd /tmp && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    cd .. &&\
    ls -ls . && \
    val3dity data/cityjson/cube.json && \
    python3 -m pytest --runfull && \
    apk del .val3dity-deps && \
    cd / && \
    apk --update add \
        make \
        gcc \
        g++ \
        boost \
        eigen \
        gmp \
        bash \
        mpfr3 && \
    rm -rf /tmp/* && \
    rm -rf /user/local/man

RUN mkdir /data && \
    chown 1001 /data && \
    chgrp 0 /data && \
    chmod g=u /data && \
    chgrp 0 /etc/passwd && \
    chmod g=u /etc/passwd
    
RUN val3dity --version

COPY --chown=1001:0 uid_entrypoint.sh /usr/local/bin/

USER 1001

WORKDIR /data

ENTRYPOINT ["/usr/local/bin/uid_entrypoint.sh"]

CMD ["--version"]
