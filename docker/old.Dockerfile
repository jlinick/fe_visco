FROM dolfinx/dolfinx
#FROM ubuntu:latest
#FROM quay.io/fenicsproject/stable:latest
MAINTAINER jlinick@mit.edu

# Set the working directory
WORKDIR /fenics
ENV HOME=${WORKDIR}
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York


RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install software-properties-common -y

# install newer python3 version
#RUN apt-get install -y \
#    wget build-essential libreadline-gplv2-dev \
#    libncursesw5-dev libssl-dev libsqlite3-dev \
#    tk-dev libgdbm-dev libc6-dev libbz2-dev \
#    libffi-dev zlib1g-dev liblzma-dev
#RUN  wget https://www.python.org/ftp/python/3.10.5/Python-3.10.5.tgz
#RUN tar xzf Python-3.10.5.tgz
#RUN cd Python-3.10.5 && ./configure --enable-optimizations
#RUN make altinstall

RUN apt-get install -y \
    zip unzip git vim wget curl

RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt update && apt install python3.9 -y
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.9 1
RUN apt remove --purge python3-apt -y
RUN apt install python3-apt -y
RUN apt install python3.9-distutils -y
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py 
RUN python3 get-pip.py --force-reinstall


RUN apt-get install -y \
    zip unzip git vim python3-pip tmux wget curl \
    software-properties-common python3-distutils \
    libgl1

#RUN pip3 install --upgrade pip && \
#    pip3 install --upgrade wheel \
#    pip3 install --upgrade setuptools \
#    pip3 install --upgrade requests \
#RUN pip3 install fenics-ffc --upgrade

RUN pip3 install \
    #numpy scipy matplotlib==3.0.3 pyproj fiona \
    numpy scipy matplotlib pyproj fiona \
    setuptools scikit-learn scikit-image tqdm pylint \
    pandas pillow plotly kaleido vedo fenics

#RUN git clone https://github.com/jtpils/vtkplotter.git
#RUN git clone https://github.com/marcomusy/vedo.git /home/fenics/vedo
#RUN cd /home/fenics/vedo && python3 setup.py install

# copy current repo
#COPY ./ /fenics


