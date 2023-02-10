FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update
RUN apt install -y python3 python3-virtualenv python3-pip
RUN pip install jupyter
RUN pip install numpy
RUN pip install matplotlib
RUN apt install -y less nano
RUN apt install -y libopenmpi3 libopenmpi-dev
RUN apt install -y libfftw3-dev libfftw3-double3 libfftw3-mpi-dev libfftw3-mpi3
RUN apt install -y build-essential
RUN apt install -y environment-modules
RUN apt install -y libnetcdf-dev libnetcdf-mpi-dev libnetcdf15 libnetcdf-mpi-13
RUN apt install -y libnetcdf-c++4-1 libnetcdf-c++4-dev
RUN apt install -y mc
ADD . /ga-aem
RUN cd /ga-aem/makefiles; ./run_make_ubuntu.sh gnu allclean 
