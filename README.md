# SSLpheno
**SSLpheno is a model for gene-phenotype association prediction by self-supervised learning**

<img src="https://github.com/bixuehua/SSLpheno/blob/main/Fig1.png">
Fig. 1.The workflow of SSLpheno. 1 Attributed network construction and feature preprocessing. 2 Self-supervised learning for pre-training. 3 Multilabel classification for phenotype prediction.

# Dependencies
Our model is implemented by Python 3.7 with Pytorch 1.4.0 and run on Nvidia GPU with CUDA 10.0

# Preprocessing
**HPO
  * Download raw gene annotation file from Human Phenotype Ontology(_http://hpo.jax.org/_), choosing the _GENES TO PHENOTYPE_ button; Click _ontology_ option under the _data_ item at the top left to download ``ontology.obo`` file.
  * Conduct ``annotation.py`` file to get raw gene set and the JSON file of gene-phenotype annotations after true_path_rule.
  * Upload the _gene_set file_ to the Uniprot website (_http://www.uniprot.org/id-mapping_) to get a mapping file from gene to protein.
**DisGeNET
  * Download ``gene_association.tsv`` file from DisGeNET (_http://www.disgenet.org/_) to get the geneSymbol.
  * Using package _disgenet2r_ provided by DisGeNET to get the gene-phenotype association in R. Genes with no disease association were filtered out. (The file is too large, we do not upload)
  * Click _UMLS CUI to several disease vocabularies_ in the _download page_ of DisGeNET to get the ``disease_mappings.tsv`` file. DisGeNET disease ids which can be mapped to HPO based on the vocabulary field in the file.
  * Based on the above steps, the gene-phenotype association file from DisGeNET was obtained.
  * ``disgenet.py`` file provide the prcessing.
**String
  * Please firstly open STRING database _(https://string-db.org/cgi/download.pl)_ and choose "organism" as "Homo sapiens", then download "9606.protein.links.v11.0.txt.gz" (version number may change). Meanwhile, download mapping file under "ACCESSORY DATA" category, or open website Uniprot website _(https://string-db.org/mapping_files/uniprot_mappings/)_ to download it. 
  * Run ``string.py`` to get a json file containing PPI data.
**Feature_GO
  * Download gene function annotation from Uniprot database _(https://string-db.org/mapping_files/uniprot_mappings/)_.
  
# src
* The implementation of SSLpheno

    ``src/evaluation.py：``This script is used to calculate macro_average and micro_average metrics.    
    ``src/layers.py：`` The module for decoding    
    ``src/model.py：`` The adaptive encoder     
    ``src/optimizer.py：`` This script is used to caculate the LOSS and optimizing model    
    ``src/main.py：``The main file of SSLpheno    
    ``src/trainNN.py：`` The deep neural network multi_label classifier    
    ``src/utils.py：``The module for loading data
	
## preprocessing
  * Python files for preprocessing raw data
  
    ``annotation.py``: The preprocessing code for _HPO_2020_.<br>
    ``disgenet.py``: The preprocessing code for _disgenet_.<br>
    ``string.py``: Providing the codes for GGA network construction.<br>
    ``obo_parser.py, ontology.py, utils.py`` are auxiliary files.
  
# Data
## hpo_2020
  * Dataset obtained from HPO database, providing the data files for SSLpheno 

    ``data/all_gene.txt：`` gene list used in our experiments   
    ``data/gene_feature_GO.txt：`` gene GO features    
    ``data/gga_string.txt：``GGAs of String    
    ``data/gene_phenotype.json：`` associations of gene-phenotype    
    ``data/gene_label.txt：``the labels of genes     
	
## disgenet
  * Dataset obtained from DisGeNET database, providing the data files for SSLpheno 
    
    ``disease_mappings.rar:``The mapping of diseases between different disease database.<br>
    ``gene_associations.tsv:``The description of genes provided in DisGeNET.<br>
    ``DisGeNet_disease_hpo.json:``The mapping of diseases and phenotypes.<br>
    ``DisGeNet_gene_hpo.json:``The association of genes and phenotypes.<br>
