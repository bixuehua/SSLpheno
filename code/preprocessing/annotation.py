#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create HPO annotations with propagation from raw file.

"""
import json
from collections import defaultdict
from ontology import HumanPhenotypeOntology
from utils import load_annotation

hpo_ontology = {
    "path": "hp_20230405.obo",
    "version": "202303"
}

config = {
  "raw_annotation": "../data/hpo_2020/genes_to_phenotype_20200825.txt",
  "ontology": "../data/hpo_2020/hp_20200811.obo",
  "processed_annotation": "../data/hpo_2020/hpo_annotation_20200825.json",
  "genes": "../data/hpo_2020/gene_list.txt"
}

# load raw hpo annotations and process
annotation = defaultdict(list)
gene_set = set()

with open(config["raw_annotation"]) as fp:
    for line in fp:
        if line.startswith('#'):
            continue
        gene_id, _, hpo_term,*_ = line.strip().split('\t')
        gene_set.add(gene_id)
        annotation[gene_id].append(hpo_term)
        annotation[gene_id] = list(set(annotation[gene_id]))

# load HPO for propagation filter in pa
ontology = HumanPhenotypeOntology(config["ontology"],version="201904")
gene_annotation = load_annotation(annotation, ontology, split=True, keep_root=True, propagate=False)
gene_annotation = gene_annotation["pa"]

# output annotation and gene
with open(config["processed_annotation"], 'w') as fp:
    json.dump(annotation, fp, indent=2)
with open(config["genes"], "w") as fp:
    for gene in gene_set:
        fp.write("%s\n" % gene)
