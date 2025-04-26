#!/bin/bash

conda conda create -n panman-env
conda activate panman-env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y panman