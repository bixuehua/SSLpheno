# SSLpheno
*SSLpheno is a model for gene-phenotype association prediction by self-supervised learning*

# Dependencies
Our model is implemented by Python 3.7 with Pytorch 1.4.0 and run on Nvidia GPU with CUDA 10.0

# Preprocessing
## HPO
  * Download raw gene annotation file from Human Phenotype Ontology(_http://hpo.jax.org/_), choosing the GENES TO PHENOTYPE button; Click _ontology_ option under the _data_ item at the top left to download _ontology.obo_ file.
  * Conduct annotation.py file to get raw gene set and the JSON file of gene-phenotype annotations after true_path_rule.
  * Upload the _gene_set file_ to the Uniprot website (_http://www.uniprot.org/id-mapping_) to get a mapping file from gene to protein.
## DisGeNET
  * Download _gene_association.tsv_ file from DisGeNET (_http://www.disgenet.org/_) to get the geneSymbol.
  * Using package _disgenet2r_ provided by DisGeNET to get the gene-phenotype association in R. Genes with no disease association were filtered out. (The file is too large, we do not upload)
  * Click _UMLS CUI to several disease vocabularies_ in the _download page_ of DisGeNET to get the _disease_mappings.tsv_ file. DisGeNET disease ids which can be mapped to HPO based on the vocabulary field in the file.
  * Based on the above steps, the gene-phenotype association file from DisGeNET was obtained.
  * _disgenet.py_ file provide the prcessing.
## String
  * Please firstly open STRING database _(https://string-db.org/cgi/download.pl)_ and choose "organism" as "Homo sapiens", then download "9606.protein.links.v11.0.txt.gz" (version number may change). Meanwhile, download mapping file under "ACCESSORY DATA" category, or open website Uniprot website _(https://string-db.org/mapping_files/uniprot_mappings/)_ to download it. 
  * Run _string.py_ to get a json file containing PPI data.
## Feature_GO
  * Download gene function annotation from Uniprot database _(https://string-db.org/mapping_files/uniprot_mappings/).
  
# src
*The implementation of SSLpheno

    src/evaluation.py：This script is used to calculate macro_average and micro_average metrics.    
    src/layers.py： The module for decoding    
    src/model.py： The adaptive encoder     
    src/optimizer.py： This script is used to caculate the LOSS and optimizing model    
    src/main.py：The main file of SSLpheno    
    src/trainNN.py： The deep neural network multi_label classifier    
    src/utils.py：The module for loading data
	
## preprocessing
  * Python files for preprocessing raw data
  
  ``annotation.py``: The preprocessing code for _HPO_2020_
  ``disgenet.py``: The preprocessing code for _disgenet_
  ``string.py``: Providing the codes for GGA network construction
  ``obo_parser.py, ontology.py, utils.py`` are auxiliary files
  
# data
## hpo_2020
  * Dataset obtained from HPO database, providing the data files for SSLpheno 

    data/all_gene_exist.txt： gene list in our experiments   
    data/edge(str)_gene_weight_GeneMANIA.txt： GGAs of GeneMANIA    
    data/edge(str)_gene_weight_HumanNet.txt：GGAs of HumanNet    
    data/edge(str)_gene_weight_String.txt： GGAs of String    
    data/hpo_annotation_20200825.json：gene-phenotype annotations       
    data/z_all_gene_(go)feature.txt：GO features of genes  
    data/z_all_gene_(hpo)label.txt：HPO labels of genes
	
## disgenet
  * Dataset obtained from DisGeNET database, providing the data files for SSLpheno 
  
	
