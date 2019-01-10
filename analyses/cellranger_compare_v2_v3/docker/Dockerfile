# Start with an ubuntu system and update it
FROM ubuntu:16.04
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list
RUN apt-get update

## Build R version 3.5.1 from source
RUN apt-get install -y build-essential
# Enable source repositories for r-base dep install
RUN apt-get build-dep -y r-base
RUN apt-get install -y wget
WORKDIR /root/
RUN wget https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz
RUN tar xzf R-3.5.1.tar.gz
WORKDIR /root/R-3.5.1
RUN ./configure --prefix=/opt/R/3.5.1/ --enable-R-shlib --with-blas --with-lapack
RUN make && make install
ENV PATH="/opt/R/3.5.1/bin/:$PATH"

# Libraries required for Seurat
RUN apt-get install -y libhdf5-dev 
RUN apt-get install -y libssl-dev

# Install git
RUN apt-get install -y git

# Install Seurat Dependencies
COPY seurat.deps.r /root/
RUN Rscript /root/seurat.deps.r

## Checkout a specific commit of Seurat
WORKDIR /root/
RUN git clone https://github.com/satijalab/seurat.git
WORKDIR /root/seurat
RUN git checkout 37a3750f599fc4c0e1226aced130441a8b250077
RUN R CMD INSTALL .

# Install pagoda2 system dependencies
RUN apt-get -y install build-essential cmake gsl-bin libgsl0-dev libeigen3-dev libboost-all-dev libssl-dev libcurl4-openssl-dev libssl-dev libcairo2-dev libxt-dev libgtk2.0-dev libcairo2-dev xvfb xauth xfonts-base

# Install pagoda2
COPY pagoda2.deps.r /root/
WORKDIR /root/
RUN Rscript pagoda2.deps.r
RUN git clone https://github.com/hms-dbmi/pagoda2.git
WORKDIR /root/pagoda2
RUN git checkout d1eadc48d5f05b1ca5f0cd6434c5d9726d7d9a8e
RUN R CMD INSTALL .

## Install gsulil
RUN apt-get install -y lsb-release curl
RUN CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg |  apt-key add -
RUN apt-get update && apt-get install -y google-cloud-sdk --allow-unauthenticated 