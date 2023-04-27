import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.modules.module import Module
from torch.nn.parameter import Parameter
import numpy as np

class SampleDecoder(Module):
    def __init__(self, act=torch.sigmoid):
        super(SampleDecoder, self).__init__()
        self.act = act

    def forward(self, zx, zy):
        sim = (zx * zy).sum(1)
        sim = self.act(sim)
    
        return sim