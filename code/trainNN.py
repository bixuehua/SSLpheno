import torch
import torch.nn as nn
from torch import optim
import torch.nn.functional as F
from torch.utils.data import DataLoader
import torch.utils.data as Data
import pandas as pd
from sklearn import metrics
import numpy as np
import sys


class Net(nn.Module):
    def __init__(self, X_train, Y_train):
        super(Net, self).__init__()
        self.dense1 = nn.Linear(X_train.shape[1], 1024)
        self.dense2 = nn.Linear(1024, 512)
        self.dense3 = nn.Linear(512, 256)
        self.dense4 = nn.Linear(256, Y_train.shape[1])
        self.dropout = nn.Dropout(0.3)
        self.LeakyReLU = nn.LeakyReLU()
        self.BatchNorm1 = nn.BatchNorm1d(1024)
        self.BatchNorm2 = nn.BatchNorm1d(512)
        self.BatchNorm3 = nn.BatchNorm1d(256)

    def forward(self, x):
        x = self.LeakyReLU(self.BatchNorm1(self.dense1(x)))
        # print(1)
        x = self.dropout(x)
        x = self.LeakyReLU(self.BatchNorm2(self.dense2(x)))
        # print(2)
        x = self.dropout(x)
        x = self.LeakyReLU(self.BatchNorm3(self.dense3(x)))
 
        x = self.dropout(x)
        out = F.softmax(F.sigmoid(self.dense4(x)), dim=1)
        return out


def train(X_train, Y_train, X_test, Y_test, lr):
    train_data = Data.TensorDataset(X_train, Y_train)
    train_loader = DataLoader(dataset=train_data, batch_size=128, shuffle=False)
    model = Net(X_train, Y_train)

    optimizer = optim.Adam(model.parameters(), lr=lr)
    for epoch in range(1000):
        s_loss = []
        for step, data in enumerate(train_loader):
            input, target = data
            # print(input.shape,target.shape)
            y_hat = model(input)
            loss = F.binary_cross_entropy(y_hat, target)
            # print('loss', loss)
            s_loss.append(loss.item())
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        sys.stdout.write('\r{} : loss = {} '.format(epoch, np.mean(s_loss)))
    y_pred = model(X_test)
    return y_pred




