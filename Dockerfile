FROM --platform=linux/amd64 mambaorg/micromamba:jammy

USER root

RUN mkdir -p /home/termite
COPY . /home/termite/
WORKDIR /home/termite
ENV PATH "/opt/conda/bin:${PATH}"
ENV PATH "/home/termite/termite_pkg:${PATH}"
ENV PYTHONPATH /home/termite

ENV HOME /home/termite


RUN apt-get update -y && apt-get install -y build-essential \
                                            libarchive-tools \
                                            curl nano git gcc \
                                            wget zlib1g \
                                            unzip zlib1g-dev
RUN micromamba install -n base -c bioconda -c anaconda -c conda-forge bedtools idr pybigwig viennarna pybedtools transtermhp samtools pysam
RUN pip install numpy==1.21 pandas scipy biopython pybigwig gffutils

RUN pip install ViennaRNA

ENTRYPOINT ["python", "termite_pkg/termite"]

