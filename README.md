# SSLpheno
SSLpheno is a model for gene-phenotype association prediction by self-supervised learning

# Dependencies
Our model is implemented by Python 3.7 with Pytorch 1.4.0 and run on Nvidia GPU with CUDA 10.0

# src
The implementation of SSLpheno
    src/evaluation.py：This script is used to calculate macro_average and micro_average metrics.
    src/layers.py： The module for decoding
    src/model.py： The adaptive encoder 
    src/optimizer.py： This script is used to caculate the LOSS and optimizing model
    src/main.py：The main file of SSLpheno
    src/trainNN.py： The deep neural network multi_label classifier
    src/utils.py：The module for loading data
    
# data
Our data files
    data/all_protein_exist.txt： Proteins from genes mapping
    data/edge(str)_protein_weight_GeneMANIA.txt： GGAs of GeneMANIA
    data/edge(str)_protein_weight_HumanNet.txt：GGAs of HumanNet
    data/edge(str)_protein_weight_String.txt： GGAs of String
    data/hpo_annotation_20200825.json：Human phenotype annotations of proteins, downloaded from https://hpo.jax.org/app/ 
    data/protein2go.json：protein-go annotations
    data/z_all_protein_(go)feature.txt：GO features of proteins
    data/z_all_protein_(hpo)label.txt：HPO labels of proteins
