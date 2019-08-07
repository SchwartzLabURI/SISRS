**Software Setup**
==================

**Software Requirements**

    The following programs must be installed and in your path prior to running SISRS:

        * Python3 (3.6.4)
            * Biopython
            * Numpy
        * Bowtie2 (2.3.4)
        * Samtools (1.3.1)
        * BBMap Suite (37.40)
        * Ray - 2.3.2-devel
        * Ray requires a working version of MPI, even if running on a single node
        * FastQC (0.11.5) (Only if trimming)

What type of machine do you plan on running *SISRS v2.0* on?

    1. Research Cluster / HPC - Clusters

        Check for the following modules on your machine:

            i. Python/3.6.4-foss-2018a
            ii. SAMtools/1..3.1-foss-2016b
            iii. bio/FastQC/0.151.5-Java-1.8.0_92
            iv. Bowtie2/2.3.4.1-intel-2018a

        Other Installs Needed:
            i. MPI - OPENMPI

               .. code-block:: bash

                   wget http://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz
                   tar xzvf openmpi-2.1.1.tar.gz
                   cd openmpi-2.1.1
                   ./configure
                   make -j4
                   make install

                   # pwd the current location
                   pwd=$(pwd)

                   # Place program into path
                   echo "export PATH=$pwd/openmpi-2.1.1:$""PATH" >> ~/.bash_profile

            ii. Ray

                .. code-block:: bash

                   git clone https://github.com/sebhtml/RayPlatform.git
                   git clone https://github.com/sebhtml/ray.git
                   cd ray
                   make PREFIX=ray-build HAVE_LIBZ=y LDFLAGS="-lpthread -lz"
                   make install

                   # pwd the current location
                   pwd=$(pwd)

                   # Place program into path
                   echo "export PATH=$pwd/ray-build:$""PATH" >> ~/.bash_profile

            iii. Biopython/Numpy

                 .. code-block:: bash

                    pip install --user biopython

                * Biopython relies heavily on numpy and auto-installs it during the biopython install

            iv. BBMap Suite

                .. code-block:: bash

                   wget -O BBMap_38.51.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_38.51.tar.gz/download \
                   && tar -xf BBMap_38.51.tar.gz \
                   && rm BBMap_38.51.tar.gz

                   # pwd the current location
                   pwd=$(pwd)

                   # build BBTools small C-lib
                   cd bbmap/jni \
                    && make -f makefile.linux

                   # Place program into path
                   echo "export PATH=$pwd/bbmap:$""PATH" >> ~/.bash_profile

            **After all of these things have been installed make sure that BBMap, Ray, and MPI is properly located in your path.**

            Entering these software packages into your path can be down through the last few lines under each code sample or it can be done through:
                .. code-block:: bash

                   vim ~/.bash_profile
                   # scroll to bottom of .bash_profile and add new line
                   Export PATH='path/to/software/package:$PATH'

        Writing A Script to Use the Modules:
            .. code-block:: bash

               #!/bin/bash
               # The script should start off will all the module and then follow with the code
               module load Python/3.6.4-foss-2018a;
               module load SAMtools/1.3.1-foss-2016b;
               module load bio/FastQC/0.11.5-Java-1.8.0_92;
               module load Bowtie2/2.3.4.1-intel-2018a

               python3 ...

    2. Running Docker or Singularity

        The Schwartz Lab @ URI is proud to also offer a docker image located at:
            docker://mcconnell22/sisrs:2.0

        **This image is an interactive image and is designed for users to run commands inside of it.**

        How to use the image:

        i. Docker

           .. code-block:: bash

              docker pull mcconnell22/sisrs:2.0
              docker run -it mcconnell22/sisrs:2.0

           This will pull the docker image and then allow you to run the docker image interactively.


        ii. Singularity

            .. code-block:: bash

               singularity shell docker://mcconnell22/sisrs:2.0

            This will pull the docker image and convert it into the appropriate singularity image. It will then launch an interactive shell immediately following.

    3. Install From Nothing

        *Refer to the section on installing software packages on research clusters to install MPI, Ray, Biopython/Numpy, and BBMap Suite*

        i. Bowtie2

           .. code-block:: bash

              wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download
              unzip download

              # Enter Bowtie2 into the PATH
              pwd=$(pwd)
              echo "export PATH=$pwd/bowtie2-2.3.3.1-linux-x86_64:$""PATH" >> ~/.bash_profile

        ii. Samtools

            .. code-block:: bash

               wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
               tar --bzip2 -xvf htslib-1.3.2.tar.bz2

               cd htslib-1.3.2
               ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
               make && \
               make install && \
               cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

               # Or any other locations
               SAMTOOLS_INSTALL_DIR=/opt/samtools

               cd /tmp
               wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
               tar --bzip2 -xf samtools-1.3.1.tar.bz2

               cd /tmp/samtools-1.3.1
               ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
               make && \
               make install

               cd /
               rm -rf /tmp/samtools-1.3.1
               ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

               # Add to the apth
               echo "export PATH=/opt/samtools/bin:$""PATH" >> ~/.bash_profile

        iii. FastQC

             .. code-block:: bash

                pwd=$(pwd)
                wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
                unzip fastqc_v0.11.5.zip
                rm fastqc_v0.11.5.zip
                chmod 755 /opt/FastQC/fastqc

                # Add to the path
                echo "export PATH=$pwd/FastQC:$""PATH" >> ~/.bash_profile
