import pickle as pkl

import networkx as nx
import numpy as np
import scipy.sparse as sp
import torch
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score
import sklearn.preprocessing as preprocess

def load_data(feature_type, net_type):
    graph = nx.Graph()

    with open("../data/all_protein_exist.txt") as f_all_protein:
           all_protein_ = f_all_protein.readlines()
        for i in range(len(all_protein_)):
            all_protein_[i] = all_protein_[i].strip("\n")
    graph.add_nodes_from(all_protein_)
    

    with open("../data/edge(str)_protein_weight_"+net_type+".txt") as f_edge:
        edge_content = f_edge.readlines()
        for line in edge_content:
            score = float(line.strip("\n").split(" ")[2])
            if score >= 0.3:
                head = line.strip("\n").split(" ")[0]  # head
                tail = line.strip("\n").split(" ")[1]  # tail
                graph.add_edge(head, tail)

    isolates_list = list(nx.isolates(graph))


    adj = nx.adjacency_matrix(graph)
    node_list = list(graph.nodes())

    if feature_type == "GO":
        all_feature_dict = {}
        with open("../data/z_all_protein_(go)feature.txt") as f_go:
            go_content = f_go.readlines()
            for line in go_content:
                protein = line.strip("\n").split(" ")[0]
                feature = line.strip("\n").split(" ")[1:]
                feature = list(map(float, feature))
                all_feature_dict[protein] = feature

        without_go_protein = []
        all_feature_dict_ = {}
        for p in node_list:
            f = all_feature_dict.get(p)
            if f == None:
                without_go_protein.append(p)
                all_feature_dict_[p] = [0.0] * 13263
            else:
                all_feature_dict_[p] = f
        print("without_go_protein:",without_go_protein)
        features = list(all_feature_dict_.values())
        features = np.array(features)

        print("PCA:")
        pca = PCA(n_components=1000)
        features = pca.fit_transform(features)

        features = torch.Tensor(features)
    if feature_type == "DC":
        all_feature_dict_3706 = {}
        with open("../data/z_all_protein_(dc)feature.txt") as f_dc:
            dc_content = f_dc.readlines()
            for line in dc_content:
                protein = line.strip("\n").split(" ")[0]
                feature = line.strip("\n").split(" ")[1:]
                feature = list(map(float, feature))
                all_feature_dict_3706[protein] = feature
        without_dc_protein = []
        all_feature_dict_ = {}
        for p in node_list:
            f = all_feature_dict.get(p)
            if f == None:
                without_dc_protein.append(p)
                all_feature_dict_[p] = [0.0] * 400
            else:
                all_feature_dict_[p] = f

        features = list(all_feature_dict_.values())
        features = np.array(features)
        features = torch.Tensor(features)

    all_hpo_dict= {}
    with open("../data/z_all_protein_(hpo)label.txt") as f_label:
        label_content = f_label.readlines()
        for line in label_content:
            protein = line.strip("\n").split(" ")[0]
            label = line.strip("\n").split(" ")[1:]
            label = list(map(int, label))
            all_hpo_dict_3709[protein] = label

    all_hpo_dict_ = {}
    for p in node_list:
        all_hpo_dict_[p] = all_hpo_dict[p]
    all_label = torch.Tensor(list(all_hpo_dict_.values()))



    return adj, features, all_label, node_list

def preprocess_graph(adj, layer, norm='sym', renorm=True):
    adj = sp.coo_matrix(adj)
    ident = sp.eye(adj.shape[0])
    if renorm:
        adj_ = adj + ident
    else:
        adj_ = adj
    
    rowsum = np.array(adj_.sum(1))
    
    if norm == 'sym':
        degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
        adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
        laplacian = ident - adj_normalized
    elif norm == 'left':
        degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -1.).flatten())
        adj_normalized = degree_mat_inv_sqrt.dot(adj_).tocoo()
        laplacian = ident - adj_normalized
        

    reg = [2/3] * (layer)

    adjs = []
    for i in range(len(reg)):
        adjs.append(ident-(reg[i] * laplacian))
    return adjs

def laplacian(adj):
    rowsum = np.array(adj.sum(1))
    degree_mat = sp.diags(rowsum.flatten())
    lap = degree_mat - adj
    return torch.FloatTensor(lap.toarray())