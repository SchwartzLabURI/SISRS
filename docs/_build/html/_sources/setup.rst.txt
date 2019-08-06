SISRS2.0 Setup
==============

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

What type of machine do you plan on running *SISRS2.0* on?

    1. Research Cluster

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

            **After all of these things have been installed make sure that BBMap, Ray, and MPI is properly located in your path. This can be done through he last few lines of each code sample.**
