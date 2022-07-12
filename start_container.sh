#!/bin/bash

# local path to this script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# name & tag of this repo to use for the image & inside the docker container
REPO='fe_visco'
TAG='0.0.1'

# output products local & docker directories
LOCAL_PRODUCT_DIR=SCRIPT_DIR #"/home/jlinick/products/${REPO}"
DOCKER_PRODUCT_DIR='/products'

# code repository path inside the docker container
DOCKER_CODE_DIR="/${REPO}"

build_dockerfile() {
    cd "${SCRIPT_DIR}"
    if [[ "$(docker images -q ${REPO}:${TAG} 2> /dev/null)" == "" ]]; then
        echo "${REPO}:${TAG} docker image does not exist, building..."
	docker build -t "${REPO}:${TAG}" -f "${SCRIPT_DIR}/docker/Dockerfile" .	
    else
        echo "${REPO}:${TAG} dockerfile exists."
    fi
}

build_dockerfile

docker run --rm -ti -v ${LOCAL_DATA_DIR}:${DOCKER_DATA_DIR}:ro -v ${SCRIPT_DIR}:${DOCKER_CODE_DIR} -v ${LOCAL_PRODUCT_DIR}:${DOCKER_PRODUCT_DIR} ${REPO}:${TAG} /${REPO}/run_simulation.py
#docker run --rm -ti --pid=host -e DISPLAY=${DISPLAY} -v /tmp/.X11-unix:/tmp/.X11-unix:ro -v ${LOCAL_DATA_DIR}:${DOCKER_DATA_DIR}:ro -v ${SCRIPT_DIR}:${DOCKER_CODE_DIR} -v ${LOCAL_PRODUCT_DIR}:${DOCKER_PRODUCT_DIR} ${REPO}:${TAG} /${REPO}/run_simulation.py

