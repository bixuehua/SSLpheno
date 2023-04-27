import numpy as np
import pandas as pd
import argparse
import os
from collections import defaultdict
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score, f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score


def get_label_frequency(ontology):
    col_sums = ontology.sum(0)
    index_11_30 = np.where((col_sums >= 11) & (col_sums <= 30))[0]
    index_31_100 = np.where((col_sums >= 31) & (col_sums<=100))[0]
    index_101_300 = np.where((col_sums>=101) & (col_sums<=300))[0]
    index_larger_300 = np.where(col_sums >= 301)[0]
    return index_11_30, index_31_100, index_101_300, index_larger_300

def f1_score_threshold(y_true, y_score, n_split=100):
    threshold = np.linspace(np.min(y_score), np.max(y_score), n_split)
    f1_list = [f1_score(y_true, y_score >= theta) for theta in threshold]
    f1_max = np.max(f1_list)
    index = f1_list.index(f1_max)
    optim_theta = threshold[index]
    return f1_max, optim_theta
    
def evaluate_performance(y_test, y_score):
    """Evaluate performance"""
    n_classes = y_test.shape[1]
    perf = dict()
    
    perf["M-aupr"] = 0.0
    perf["M-fmax"] = 0.0
    perf["M-auc"] = 0.0
    n = 0
    aupr_list = []
    num_pos_list = []
    y_score_all = np.zeros(0)
    y_test_all = np.zeros(0)
    for i in range(n_classes):
        num_pos = sum(y_test[:, i])
        if len(np.unique(y_test[:, i])) < 2:
            continue
        # if num_pos > 1: ### BXH
        if num_pos > 0: ### BXH
            auc = compute_roc(y_test[:, i], y_score[:, i])
            perf["M-auc"] += auc
            ap = average_precision_score(y_test[:, i], y_score[:, i])
            n += 1
            perf["M-aupr"] += ap
            aupr_list.append(ap)
            num_pos_list.append(num_pos)
            
            fmax, _ = f1_score_threshold(y_test[:, i], y_score[:, i])
            perf["M-fmax"] += fmax
            y_score_all = np.concatenate((y_score_all, y_score[:, i]))
            y_test_all = np.concatenate((y_test_all, y_test[:, i]))
    print('n', n,len(y_test_all))
    perf["M-aupr"] /= n
    perf["M-auc"] /= n
    perf["M-fmax"] /= n

    # Compute micro-averaged AUPR
    perf['m-aupr'] = average_precision_score(y_test_all, y_score_all)
    perf['m-fmax-hpodnets'], perf['optim_theta'] = f1_score_threshold(y_test_all, y_score_all)

    return perf

def compute_roc(labels, preds):
    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    roc_auc = auc(fpr, tpr)
    return roc_auc

def get_results(ontology, Y_test, y_score):
    perf = defaultdict(dict)

    index_11_30, index_31_100, index_101_300, index_301 = get_label_frequency(ontology)
    print('index 11-30', len(index_11_30),len(index_31_100),len(index_101_300),len(index_301))
    perf['11-30'] = evaluate_performance(Y_test[:,index_11_30], y_score[:,index_11_30])
    perf['31-100'] = evaluate_performance(Y_test[:,index_31_100], y_score[:,index_31_100])
    perf['101-300'] = evaluate_performance(Y_test[:,index_101_300], y_score[:,index_101_300])
    perf['301-'] = evaluate_performance(Y_test[:,index_301], y_score[:,index_301])
    perf['all'] = evaluate_performance(Y_test, y_score)
    perf = pd.DataFrame(perf)
    
    return perf


if __name__ == "__main__":

    ### load result y_score
    y_score = pd.read_pickle('result/after_precessfeature_test_pred_mlp_500_500.pkl')

    ### load test label
    graph, num_features, num_classes, node_list = load_my_data()
    test_mask = graph.ndata["test_mask"]
    labels = graph.ndata["label"]

    y_test = labels[test_mask]

    y_test = y_test.cpu().detach().numpy().astype(float)
    perf = get_results(labels, y_test, y_score)
    print(type(perf))


    with open("result.csv", "w") as f:
        perf.to_csv(f)