FROM pymor/fenics_py3.9
MAINTAINER jlinick@mit.edu

# Set the working directory
WORKDIR /fe_visco
ENV HOME=${WORKDIR}
RUN apt-get update -y && apt-get upgrade -y

RUN apt-get install -y \
    zip unzip git vim python3-pip tmux wget curl \
    software-properties-common python3-distutils \
    libgl1

RUN pip3 install \
    numpy scipy matplotlib pyproj fiona \
    setuptools scikit-learn scikit-image tqdm pylint
#pandas pillow plotly kaleido vedo fenics

ENTRYPOINT ["/bin/sh", "-c"]
CMD ["bash"]
