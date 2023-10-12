#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct STRING PPI network.

Open website https://string-db.org/cgi/download.pl and choose organism as
Homo sapiens, then download interaction data "9606.protein.links.v11.0.txt.gz".
Besides, download mapping file under "ACCESSORY DATA" category, or open website
https://string-db.org/mapping_files/uniprot_mappings/ to download it.
"""
import json
from collections import defaultdict
from utils import gene2uniprot

def string2uniprot(file_path):
    """Map STRING accession to UniProt accession.
    :param file_path: path to mapping file provided by STRING website
        url: https://string-db.org/mapping_files/uniprot_mappings/
    :return: dict, key: STRING accession, value: UniProt accession
        { string_ac1: uniprot_ac1, string_ac2: uniprot_ac2, ... }
    """
    mapping = dict()
    with open(file_path) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            entries = line.strip().split('\t')
            uniprot_ac = entries[1].split('|')[0]
            string_ac = entries[2] if entries[2].startswith('9606.') \
                else '9606.' + entries[2]
            mapping[string_ac] = uniprot_ac
    return mapping


def get_string_network(path_to_network, path_to_mapping):
    """Construct STRING PPI network.
    :param path_to_network: path to STRING network data
        url: https://string-db.org/cgi/download.pl
    :param path_to_mapping: path to mapping file
        url: https://string-db.org/mapping_files/uniprot_mappings/
    :return: dict, PPI network
        { protein1: { protein1a: score1a, protein1b: score1b, ... },
          protein2: { protein2a: score2a, protein2b: score2b, ... },
          ... }
    """
    network = defaultdict(dict)
    mapping = string2uniprot(path_to_mapping)
    with open(path_to_network) as fp:
        for line in fp:
            if line.startswith("protein1"):
                continue
            string_ac1, string_ac2, score = line.strip().split()
            try:    # if no matched accession found, pass it
                protein1 = mapping[string_ac1]
                protein2 = mapping[string_ac2]
            except KeyError:
                continue
            score = float(score) / 1000
            network[protein1][protein2] = network[protein2][protein1] = score
    return network


if __name__ == "__main__":

    config = {
        "gene_list": "../hpo_2020/gene_sets.txt",
        "network": "../network/9606.protein.links.v11.0.txt",
        "mapping": "../network/human.uniprot_2_string.2018.tsv",
        "gene_to_protein": "../network/gene2uniprot_20200825.txt",
        "ppi": "../network/STRING.v11.0.json",
        "gga": "../network/gga.txt"
    }

    # # load mapping of gene id to uniprot id
    protein2gene = gene2uniprot(config["gene_to_protein"], gene_column=1, uniprot_column=0)

    # get GGA network
    protein_mapping = string2uniprot(config["mapping"])
    network = defaultdict(dict)
    gene2gene = defaultdict(dict)
    proteins = set(protein2gene.keys())
    with open(config['network']) as fp:
        for line in fp:
            if line.startswith("protein1"):
                continue
            string_ac1, string_ac2, score = line.strip().split()
            try:  # if no matched accession found, pass it
                protein1 = protein_mapping[string_ac1]
                protein2 = protein_mapping[string_ac2]
            except KeyError:
                continue
            score = float(score) / 1000
            if protein1 in proteins and protein2 in proteins:
                for each in protein2gene[protein1]:
                    for item in protein2gene[protein2]:
                        gene2gene[each][item] = gene2gene[each][item] = score

    # write into file
    with open(config["gga"], 'w') as fp:
        for each in gene2gene.keys():
            for item in gene2gene[each].keys():
                fp.write(each + ' ' + item + ' ' + str(gene2gene[each][item]) + '\n')
