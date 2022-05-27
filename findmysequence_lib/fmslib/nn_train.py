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

from optparse import OptionParser, SUPPRESS_HELP

import torch
import torch.nn.functional as F
from torch import optim
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset

import numpy as np

from fmslib import nn_utils
from fmslib import xyz_utils
from fmslib import nn_models



def parse_args():
    """setup program options parsing"""
    parser = OptionParser()


    parser.add_option("--data", action="store", \
                            dest="data", type="string", metavar="FILENAME", \
                  help="Input HDF5(s) with training data", default=None)

    parser.add_option("-o", action="store", dest="ofname", type="string", metavar="FILENAME", \
                  help="output model", default=None)

    parser.add_option("--epochs", action="store", dest="epochs", type="int", metavar="INT", \
                  help="training epochs [default: %default]", default=100)

    parser.add_option("--nodes", action="store", dest="nodes", type="int", metavar="INT", \
                  help="mid-layer nodes [default: %default]", default=None)


    parser.add_option("--lr", action="store", dest="lr", type="float", metavar="FLOAT", \
                  help="learning rate [default: %default]", default=1e-4)

    parser.add_option("--valfrac", action="store", dest="valfrac", type="float", metavar="FLOAT", \
                  help="calidation set fraction [default: %default]", default=0.1)

    parser.add_option("--validate", action="store", dest="validate", type="int", metavar="INT", \
                  help="number of validation cycles [default: %default]", default=0)

    parser.add_option("--test", action="store_true", dest="test", default=False, \
                  help="a simple test")

    parser.add_option("--skladaczka", action="store_true", dest="skladaczka", default=False, \
                  help=SUPPRESS_HELP)

    (options, _args)  = parser.parse_args()
    return (parser, options)


def accuracy(out, yb):
    preds = torch.argmax(out, dim=1)
    return (preds == yb).float().mean()



def train(d_train, l_train, d_valid, l_valid, epochs=100, lr=1e-3, model_ofname=None, nodes=None):

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    print("device: ", device)


    xyzu = xyz_utils.xyz_utils()


    loss_func = F.cross_entropy


    x_train = torch.tensor(d_train, device=device).float()
    y_train = torch.tensor(l_train, device=device).long()

    x_valid = torch.tensor(d_valid, device=device).float()
    y_valid = torch.tensor(l_valid, device=device).long()


    if nodes is None:
        print("Using naive model")
        model = nn_models.naive().to(device)
    else:
        print("Using model with %i hidden nodes"%nodes)
        model = nn_models.simple(nn=nodes).to(device)

    model_params = 0
    for parameter in model.parameters():
        _p=1
        for _n in parameter.size(): _p*=_n
        model_params+=_p
    print("Total number of model parameters: %i" % model_params)
    print()




    opt = optim.SGD(model.parameters(), lr=0.01)
    opt = optim.Adam(model.parameters(), lr=lr, weight_decay=0.1*lr)


    bs = 200
    train_ds = TensorDataset(x_train, y_train)
    train_dl = DataLoader(train_ds, batch_size=bs)

    valid_ds = TensorDataset(x_valid, y_valid)
    valid_dl = DataLoader(valid_ds, batch_size=bs * 2)

    for epoch in range(epochs):
        model.train()
        for _xb, _yb in train_dl:

            xb = _xb.to(device)
            yb = _yb.to(device)

            pred = model(xb)
            loss = loss_func(pred, yb)

            loss.backward()
            opt.step()
            opt.zero_grad()

        model.eval()
        with torch.no_grad():
            valid_loss = sum(loss_func(model(xb), yb) for xb, yb in valid_dl)

        #print(epoch, valid_loss / len(valid_dl))
        print("%-4i"%epoch, "pen_valid= %f"%loss_func(model(x_valid), y_valid),
                            "pen_train= %f"%loss_func(model(x_train), y_train),
                            "acc_train= %f"%accuracy(model(x_train), y_train),
                            "acc_valid= %f"%accuracy(model(x_valid), y_valid))

        if epoch %1000 == 0:
            torch.save(model.state_dict(), 'tmp_model.dat' if model_ofname is None else model_ofname)

            for ii in range(20):
                print( "%3s  %5.2f" % (xyzu.standard_aa_names[ii],
                                          accuracy(model(x_train[y_train==ii]), y_train[y_train==ii]).float()) )




    #model.load_state_dict(torch.load('tmp_model.dat'))
    #model.eval()
    #print(epoch, loss_func(model(x_valid), y_valid), accuracy(model(x_train), y_train), accuracy(model(x_valid), y_valid))
    d = {}
    for ii in range(20):
        d[xyzu.standard_aa_names[ii]] = accuracy(model(x_train[y_train==ii]), y_train[y_train==ii]).item()

    return d


def complete_training(input_datafiles, model_ofname=None, epochs=500, validate=0, lr=1e-5, valfrac=0.1, nodes=None):

    nnu = nn_utils.nn_utils()
    idata, ilabels = nnu.read_multiple_sets(input_datafiles)
    data = {}

    for ii in range(validate+1):
        d_train, l_train, d_valid, l_valid = nnu.select_trainset(idata, ilabels, valid_fraction=valfrac)
        print("Selected %i/training and %i/validation" % (len(l_train), len(l_valid)))
        d = train(d_train, l_train, d_valid, l_valid, epochs=epochs, lr=lr, model_ofname=model_ofname, nodes=nodes)
        for k,v in d.items():
            _a = data.setdefault(k, [])
            _a.append(v)

    for k,v in data.items():
        print( "%3s%10.3f%10.3f" % (k, np.mean(v), np.std(v)) )


def basic_test():
    """
        It should result in a nice, indetectable overfitting
    """

    complete_training(["tst_prot_data_6dg7_refined.hdf5"], epochs=100, lr=1e-3)


def skladaczka():

    import torch
    import nn_models

    model = nn_models.naive()
    fname=os.path.join('..', 'data', 'nn_gpu_1e4_bs200_naive.dat')
    #fname=os.path.join('..', 'data', 'xtst_5k.dat')
    model.load_state_dict(torch.load(fname, map_location=torch.device('cpu')))
    model=model#.double()
    model.eval()
    example = torch.zeros(1, 324)
    module = torch.jit.trace(model, example)#.double()
    for i in range(1):
        inp = np.random.rand(324)
        inp-=np.mean(inp)
        inp/=np.std(inp)

        output = module(torch.tensor(torch.ones(1, 324)))
        print(output)
    sm = torch.jit.script(module)
    module.save("em_model.pt")

def main():

    (parser, options) = parse_args()

    print( " ==> Command line: %s" % (" ".join(sys.argv)) )

    if options.test:
        basic_test()
        exit(0)

    if options.skladaczka:
        skladaczka()
        exit(0)

    if options.data is None:
        parser.print_help()
        exit(1)


    # --data must accept wildcards and will be processed separately: create list with input hdf5 file names
    input_datafiles = []
    _idx = sys.argv.index("--data")+1

    while _idx<len(sys.argv) and not sys.argv[_idx].startswith("-"):
        input_datafiles.append( sys.argv[_idx] )
        _idx+=1

    complete_training(input_datafiles,
                      epochs=options.epochs,
                      validate=options.validate,
                      lr=options.lr,
                      valfrac=options.valfrac,
                      model_ofname=options.ofname,
                      nodes=options.nodes)


if __name__=="__main__":
    main()
