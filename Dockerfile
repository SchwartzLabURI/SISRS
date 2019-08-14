# ubuntu OS --> from Scratch
FROM ubuntu

SHELL ["/bin/bash", "-c"]
RUN apt-get update && apt-get install -y python3 \
    python \
    python3-pip \
    python-pip \
    apt-utils \
    bzip2 \
    gcc \
    make \
    ncurses-dev \
    wget \
    zlib1g-dev \
    htop \
    locate \
    unzip \
    git \
    vim \
    screen \
    libpthread-stubs0-dev \
    perl \
    locales

RUN pip3 install biopython && pip3 install pandas && pip install biopython

# Ad python to the path

# Set the locale
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME


# Setup bowtie2
RUN mkdir /opt/bowtie2/
WORKDIR /opt/bowtie2/

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download
RUN unzip download

# Setup Samtools
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2

WORKDIR /tmp/htslib-1.3.2
RUN ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar --bzip2 -xf samtools-1.3.1.tar.bz2

WORKDIR /tmp/samtools-1.3.1
RUN ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install

WORKDIR /
RUN rm -rf /tmp/samtools-1.3.1
RUN ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

# setup OPENMPI
WORKDIR /opt
RUN wget http://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz
RUN tar xzvf openmpi-2.1.1.tar.gz
WORKDIR /opt/openmpi-2.1.1
RUN ./configure
RUN make -j4
RUN make install
RUN ldconfig
RUN rm /opt/openmpi-2.1.1.tar.gz

# setup RAY
WORKDIR /opt
RUN git clone https://github.com/sebhtml/RayPlatform.git
RUN git clone https://github.com/sebhtml/ray.git
WORKDIR /opt/ray
RUN make PREFIX=ray-build HAVE_LIBZ=y LDFLAGS="-lpthread -lz"
RUN make install

# setup FASTQC
WORKDIR /opt
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
RUN unzip fastqc_v0.11.5.zip
RUN rm fastqc_v0.11.5.zip
RUN chmod 755 /opt/FastQC/fastqc

# setup BBMAP
WORKDIR /opt
RUN BBMAP_VERSION=38.51 \
    && BBMAP=BBMap_$BBMAP_VERSION.tar.gz \
    && wget -O $BBMAP https://sourceforge.net/projects/bbmap/files/$BBMAP/download \
    && tar -xf $BBMAP \
    && rm $BBMAP


# build BBTools small C-lib
RUN cd /opt/bbmap/jni \
    && make -f makefile.linux

# Add everything to the path
ENV PATH "/opt/bowtie2/bowtie2-2.3.3.1-linux-x86_64:/opt/samtools/bin:/opt/openmpi-2.1.1:/opt/ray/ray-build:/opt/FastQC:/opt/bbmap:$PATH"

RUN adduser --home /home/sisrs --disabled-password --disabled-login sisrs

USER sisrs

WORKDIR /home/sisrs

ENTRYPOINT ["/bin/bash"]
