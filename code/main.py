from __future__ import division
from __future__ import print_function
import os, sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
# For replicating the experiments
SEED = 42
import argparse
import time
import random
import numpy as np
import scipy.sparse as sp
import torch
import json

np.random.seed(SEED)
torch.manual_seed(SEED)
from torch import optim
import torch.nn.functional as F
from model import LinTrans, LogReg
from optimizer import loss_function
from utils import *
from tqdm import tqdm
from sklearn.preprocessing import normalize


from trainNN import train_nn
from evaluation_macro import get_results

parser = argparse.ArgumentParser()
parser.add_argument('--gnnlayers', type=int, default=4, help="Number of gnn layers")
parser.add_argument('--linlayers', type=int, default=1, help="Number of hidden layers")
parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')
parser.add_argument('--dims', type=int, default=[2048], help='Number of units in hidden layer 1.')
parser.add_argument('--lr', type=float, default=0.001, help='Initial learning rate.')
parser.add_argument('--upth_st', type=float, default=0.011, help='Upper Threshold start.')
parser.add_argument('--lowth_st', type=float, default=0.1, help='Lower Threshold start.')
parser.add_argument('--upth_ed', type=float, default=0.001, help='Upper Threshold end.')
parser.add_argument('--lowth_ed', type=float, default=0.5, help='Lower Threshold end.')
parser.add_argument('--upd', type=int, default=10, help='Update epoch.')
parser.add_argument('--bs', type=int, default=10000, help='Batchsize.')
parser.add_argument('--dataset', type=str, default='citeseer', help='type of dataset.')
parser.add_argument('--no-cuda', action='store_true', default=False, help='Disables CUDA training.')
args = parser.parse_args()
args.cuda = not args.no_cuda and torch.cuda.is_available()
print("args.cuda：",args.cuda)

if args.cuda is True:
    print('Using GPU')
    torch.cuda.manual_seed(SEED)
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"


def update_similarity(z, upper_threshold, lower_treshold, pos_num, neg_num):
    f_adj = np.matmul(z, np.transpose(z))
    cosine = f_adj
    cosine = cosine.reshape([-1,])
    pos_num = round(upper_threshold * len(cosine))
    neg_num = round((1-lower_treshold) * len(cosine))
    
    pos_inds = np.argpartition(-cosine, pos_num)[:pos_num]
    neg_inds = np.argpartition(cosine, neg_num)[:neg_num]
    
    return np.array(pos_inds), np.array(neg_inds)

def update_threshold(upper_threshold, lower_treshold, up_eta, low_eta):
    upth = upper_threshold + up_eta
    lowth = lower_treshold + low_eta
    return upth, lowth


def gae_for(args):
    # net_list = ["String", "HumanNet", "Gene"]
    net_list = ["String"]
    all_mu = []
    for net in net_list:
        adj, features, true_labels, node_list = load_data(feature_type="GO", net_type=net)

        n_nodes, feat_dim = features.shape
        dims = [feat_dim] + args.dims
    
        layers = args.linlayers
    
        adj = adj - sp.dia_matrix((adj.diagonal()[np.newaxis, :], [0]), shape=adj.shape)
        adj.eliminate_zeros()

    
        n = adj.shape[0]
    
        adj_norm_s = preprocess_graph(adj, args.gnnlayers, norm='sym', renorm=True)
        sm_fea_s = sp.csr_matrix(features).toarray()
    
        print("The type of GGA："+net+'，Laplacian Smoothing...')
        for a in adj_norm_s:
            sm_fea_s = a.dot(sm_fea_s)



        adj_1st = (adj + sp.eye(n)).toarray()
    
    
        adj_label = torch.FloatTensor(adj_1st)
    
        model = LinTrans(layers, dims)
    
        optimizer = optim.Adam(model.parameters(), lr=args.lr)
    
        sm_fea_s = torch.FloatTensor(sm_fea_s)
        adj_label = adj_label.reshape([-1,])
    
        if args.cuda:
            model.cuda()
            inx = sm_fea_s.cuda()
            adj_label = adj_label.cuda()
        else:
            inx = sm_fea_s
            adj_label = adj_label
    
        pos_num = len(adj.indices)
        neg_num = n_nodes*n_nodes-pos_num

        up_eta = (args.upth_ed - args.upth_st) / (args.epochs/args.upd)
        low_eta = (args.lowth_ed - args.lowth_st) / (args.epochs/args.upd)

        pos_inds, neg_inds = update_similarity(normalize(sm_fea_s.numpy()), args.upth_st, args.lowth_st, pos_num, neg_num)
        upth, lowth = update_threshold(args.upth_st, args.lowth_st, up_eta, low_eta)

        bs = min(args.bs, len(pos_inds))
        length = len(pos_inds)

        pos_inds_cuda = torch.LongTensor(pos_inds).cuda()
        print('Start Training...')
        for epoch in tqdm(range(args.epochs)):

            st, ed = 0, bs
            batch_num = 0
            model.train()
            length = len(pos_inds)

            while ( ed <= length ):
                sampled_neg = torch.LongTensor(np.random.choice(neg_inds, size=ed-st)).cuda()
                sampled_inds = torch.cat((pos_inds_cuda[st:ed], sampled_neg), 0)
                t = time.time()
                optimizer.zero_grad()
                xind = sampled_inds // n_nodes
                yind = sampled_inds % n_nodes
                x = torch.index_select(inx, 0, xind)
                y = torch.index_select(inx, 0, yind)
                zx = model(x)
                zy = model(y)
                batch_label = torch.cat((torch.ones(ed-st), torch.zeros(ed-st))).cuda()
                batch_pred = model.dcs(zx, zy)
                loss = loss_function(adj_preds=batch_pred, adj_labels=batch_label, n_nodes=ed-st)

                loss.backward()
                cur_loss = loss.item()
                optimizer.step()

                st = ed
                batch_num += 1
                if ed < length and ed + bs >= length:
                    ed += length - ed
                else:
                    ed += bs

            if (epoch + 1) % args.upd == 0:
                model.eval()
                mu = model(inx)
                hidden_emb = mu.cpu().data.numpy()
                upth, lowth = update_threshold(upth, lowth, up_eta, low_eta)
                pos_inds, neg_inds = update_similarity(hidden_emb, upth, lowth, pos_num, neg_num)
        mu = model(inx)
        all_mu.append(mu)

    if len(net_list) == 1:
        final_mu = all_mu[0]
    if len(net_list) == 3:
        final = []
        for i in range(3709):
            final.append(all_mu[0][i].cpu().detach().numpy().tolist()+all_mu[1][i].cpu().detach().numpy().tolist()+all_mu[2][i].cpu().detach().numpy().tolist())
        final_mu = torch.tensor(final)
        print(final_mu.shape)
    
    
    
    all_X = final_mu
    all_Y = true_labels

    all_X = all_X.cpu().detach().numpy()
    all_Y = all_Y.cpu().detach().numpy()
    
    print(all_X.shape, type(all_X), all_Y.shape, type(all_Y))

    num_test = int(np.floor(all_Y.shape[0] / 10.))
    num_train = all_Y.shape[0] - num_test
    all_idx = list(range(all_Y.shape[0]))
    
    np.random.shuffle(all_idx)

    train_idx = all_idx[:num_train]
    test_idx = all_idx[num_train:(num_train + num_test)]

    X_train = all_X[train_idx]
    X_test = all_X[test_idx]

    Y_train_hpo = all_Y[train_idx]
    Y_test_hpo = all_Y[test_idx]

    print("Train DNN classifier！！！！")
    print(X_train.shape,Y_train_hpo.shape,X_test.shape,Y_test_hpo.shape)
    y_score_hpo = train_nn(X_train, Y_train_hpo, X_test, Y_test_hpo)
    

    print("y_score_hpo.shape：", y_score_hpo.shape)
    print(type(all_Y), type(Y_test_hpo), type(y_score_hpo))
    print(all_Y.shape)
    print(Y_test_hpo.shape)
    print(y_score_hpo.shape)
    perf_hpo = get_results(all_Y, Y_test_hpo, y_score_hpo)
    print(perf_hpo)


if __name__ == '__main__':
    gae_for(args)

