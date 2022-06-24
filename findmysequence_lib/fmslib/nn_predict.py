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
__date__ = "14 Dec 2021"


import os,sys,re

from optparse import OptionParser, SUPPRESS_HELP

ROOT = os.path.dirname(os.path.realpath(__file__))

import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np

USE_TORCH = 0

try:
    _ = ctypes.CDLL(os.path.join(ROOT, 'fms_nnmodel.so'))
except:
    import torch
    import torch.nn.functional as F
    from torch import optim
    from torch.utils.data import DataLoader
    from torch.utils.data import TensorDataset

    from fmslib import nn_models
    USE_TORCH = 1

    print(" ==> Loaded native pytorch module")

import numpy as np

#from fmslib import nn_utils
#from fmslib import xyz_utils
from fmslib import descriptors


class nn_predict:


    def __init__(self, mapin=None, mtzin=None, labin=None, badd=None):
        if mapin is None:
            self.__load_model(em=False)
        else:
            self.__load_model(em=True)

        self.descr = descriptors.descriptors(mapin=mapin, mtzin=mtzin, labin=labin, badd=badd)


    def __load_model(self, em=True):
        self.em_model = em

        if self.em_model:
            print(" ==> Parsed EM model")
            fname=os.path.join(ROOT, '..', 'data', 'nn_gpu_1e4_bs200_naive.dat')
            tscript_fname=os.path.join(ROOT, 'em_model.pt')
        else:
            print(" ==> Parsed XRAY model")
            fname=os.path.join(ROOT, '..', 'data', 'xtst_5k.dat')
            tscript_fname=os.path.join(ROOT, 'xray_model.pt')

        if USE_TORCH:
            device = torch.device("cpu")

            self.model = nn_models.naive().to(device)

            self.model.load_state_dict(torch.load(fname, map_location=device))
            self.model.eval()
        else:
            self.fms_nnlib = ctypes.CDLL(os.path.join(ROOT, 'fms_nnmodel.so'))
            self.fms_nnlib.predict.restype = ndpointer(dtype=ctypes.c_float, shape=(20,))


    #def tst(self):

    #    nnu = nn_utils.nn_utils()
    #    idata, ilabels = nnu.read_multiple_sets([os.path.join(ROOT, "tst_prot_data_6dg7_refined.hdf5")])

    #    pred = torch.argmax( self.model(torch.tensor(idata).float()), dim=1)
    #    true = torch.tensor(ilabels).long()
    #    print( (pred == true) .float().mean() )
    #    nn,tr=0,0
    #    for l,d in zip(ilabels, idata):
    #        if l == torch.argmax(self.model(torch.tensor([d]).float()), dim=1 ): tr+=1
    #        nn+=1
    #    print(nn, tr/nn)


    #def tst_predict(self, modelin):

    #    ph, symm = self.descr._xyz_utils.read_ph(modelin)
    #    nn=0
    #    tr=0

    #    for res in ph.atom_groups():
    #        res_type, cloud = self.descr.describe(res)
    #        if res_type is None: continue

    #        idx = torch.argmax( self.model(torch.tensor([cloud]).float()), dim=1)
    #        if self.descr._xyz_utils.standard_aa_names.index(res_type)==idx: tr+=1
    #        nn+=1
    #    print(nn, tr/nn)




    def predict(self, res):

        res_type, cloud = self.descr.describe(res)
        if res_type is None: return None
        if USE_TORCH:
            pred = self.model(torch.tensor(np.array([cloud])).float()).exp().detach().numpy().flatten()
        else:
            #pred = np.exp(self.model.predict(cloud.astype('f'))[0])
            FloatArray324 = ctypes.c_float * 324
            parameters_array = FloatArray324(*cloud)
            pred = self.fms_nnlib.predict(parameters_array, 1 if self.em_model else 0)

        pred_dict = {}

        for r, p in zip(self.descr._xyz_utils.standard_aa_names, pred):
            pred_dict[r] = p
        return pred_dict




def test_nn(n=10):
    import torch

    fname=os.path.join(ROOT, '..', 'data', 'nn_gpu_1e4_bs200_naive.dat')
    device = torch.device("cpu")

    from fmslib import nn_models
    pt_model = nn_models.naive().to(device)
    pt_model.load_state_dict(torch.load(fname, map_location=device))
    pt_model.eval()

    ct_model = ctypes.CDLL(os.path.join(ROOT, 'fms_nnmodel.so'))
    ct_model.predict.restype = ndpointer(dtype=ctypes.c_float, shape=(20,))

    for ii in range(n):
        cloud = np.random.normal(10, 10.0, 324)
        cloud = (cloud-cloud.mean())/cloud.std()


        pt_pred = pt_model(torch.tensor(np.array([cloud])).float()).exp().detach().numpy().flatten()

        FloatArray324 = ctypes.c_float * 324
        _array = FloatArray324(*cloud)
        ct_pred = ct_model.predict(_array, 1)

        err = np.linalg.norm(ct_pred-pt_pred)
        assert err<1e-3, (pt_pred, ct_pred)


if __name__=="__main__":
    test_nn(n=10000)
