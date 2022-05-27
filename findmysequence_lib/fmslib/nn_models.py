#! /usr/bin/env libtbx.python

# Copyright (c) 2021, EMBL/Grzegorz Chojnowski (gchojnowski@embl-hamburg.de)
# All rights reserved.
# License: BSD 3-Clause License

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

__author__ = "Grzegorz Chojnowski"
__date__ = "16 Feb 2021"

import os,sys,re

import torch
import torch.nn.functional as F
from torch import optim
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset




class naive(torch.nn.Module):

    def __init__(self):
        super().__init__()

        self.dense1 = torch.nn.Linear(in_features=324, out_features=324)
        self.dense2 = torch.nn.Linear(in_features=324, out_features=20)
        self.do = torch.nn.Dropout(p=0.5)

    def forward(self, state):

        #return F.log_softmax(self.dense2(state), dim=1)

        x = F.relu(self.dense1(state))
        x1 = self.do(x)
        outputs = self.dense2(x1)
        return F.log_softmax(outputs, dim=1)

class simple(torch.nn.Module):

    def __init__(self, nn=100):
        super().__init__()

        self.dense1 = torch.nn.Linear(in_features=324, out_features=nn)
        self.dense2 = torch.nn.Linear(in_features=nn, out_features=20)
        self.do = torch.nn.Dropout(p=0.5)

    def forward(self, state):

        x = F.relu(self.dense1(state))
        x1 = self.do(x)
        outputs = self.dense2(x1)
        return F.log_softmax(outputs, dim=1)

class simple_2l(torch.nn.Module):

    def __init__(self, nn=100):
        super().__init__()

        self.dense1 = torch.nn.Linear(in_features=324, out_features=nn)
        self.dense2 = torch.nn.Linear(in_features=nn, out_features=100)
        self.dense3 = torch.nn.Linear(in_features=100, out_features=20)
        self.do1 = torch.nn.Dropout(p=0.5)
        self.do2 = torch.nn.Dropout(p=0.5)

    def forward(self, state):

        x = F.relu(self.dense1(state))
        x1 = self.do1(x)
        x2 = F.relu(self.dense2(x1))
        x3 = self.do2(x2)
        outputs = self.dense3(x3)

        return F.log_softmax(outputs, dim=1)

