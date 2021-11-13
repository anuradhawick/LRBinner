FROM continuumio/miniconda3
LABEL maintainer="anruadhawick@gmail.com"

RUN conda install -y pytorch cudatoolkit=11.3 -c pytorch
RUN conda install -c conda-forge biopython
RUN conda install -y tqdm tabulate
RUN conda install -y -c anaconda scikit-learn
RUN apt-get update && \
    apt-get -y install g++ bzip2 lzma-dev zlib1g-dev libbz2-dev && \
    apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev liblzma-dev libomp-dev

RUN mkdir /usr/LRBinner
COPY . /usr/LRBinner/
WORKDIR /usr/LRBinner/

RUN ["bash", "docker_build.sh"]

ENV PATH="/usr/LRBinner/:${PATH}"

ENTRYPOINT ["LRBinner"]
CMD ["--help"]

# docker build -t anuradhawick/lrbinner .
# docker run  --rm -it --gpus all -v `pwd`:`pwd` -u `id -u`:`id -g`  anuradhawick/lrbinner
# docker run  --rm -it --gpus all -v `pwd`:`pwd` -u `id -u`:`id -g` --entrypoint /bin/bash  anuradhawick/lrbinner
