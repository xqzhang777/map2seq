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
__date__ = "22 Oct 2021"

import os,sys,re

input_em="""
-0.2681321,1.0402813,-0.019997548,0.8053067,0.26752535,-1.9855213,-0.7540464,-0.16831578,-1.4197662,-1.9409006,-1.7499865,-0.57322794,0.20964006,-0.91138095,-0.7110432,-0.59838814,-0.021187957,-0.29668984,-0.8859751,0.11075916,-0.41084698,-1.0017552,-0.42335895,-0.5094864,-0.2759891,-2.124697,-1.8803382,-1.2751404,-0.37930968,-1.643896,-1.468796,-1.4691317,-0.70821667,0.9117377,1.3653181,0.93195546,0.21707126,-0.34669924,0.11062859,-0.778974,-0.01577638,-0.7008413,-1.4440293,-1.1275023,-0.678824,-0.27194768,0.13840172,-0.4958808,-0.5910657,2.6167946,1.7783097,0.4665805,-0.054516103,-0.47718555,1.6534646,1.5617867,1.3627502,0.90162903,0.4027058,-0.7253437,0.087974034,0.2729081,0.49885142,0.764527,0.64844745,-1.2157708,-1.4388376,-0.9397317,-0.51515174,0.046122998,0.1429022,1.1521883,0.8582204,0.09746878,-0.75084627,-0.47969607,-0.34859952,-0.6572948,2.4160738,1.196757,0.1492668,2.3040617,1.0371252,-0.34144258,-0.6278315,-0.042191345,-0.17378883,-1.2326062,-0.4094434,-0.016610237,-1.3552305,-0.88448673,-0.2660648,-0.30057144,-0.05011212,0.0063696736,-0.017802553,-0.05249696,-0.401318,3.0850465,2.300961,2.0666842,1.3227632,0.18296519,-0.45068017,2.9080343,3.661098,2.9770555,2.2392538,1.2791339,0.34440127,1.04705,1.7420827,1.694114,1.5290825,1.4256687,0.576179,-0.38574892,-0.66770023,-0.70405346,0.2007557,0.53811896,0.28137755,2.0912926,0.38103446,-1.1121341,-1.6880342,-0.99096245,-0.039083626,0.22595356,2.7812278,0.121447176,-1.2455597,-0.9563433,0.43749687,-0.9376573,0.5062168,0.55654526,-0.3932322,-0.1910638,0.33372188,-0.9949559,-0.36572576,0.010702438,-0.19771121,-0.10163034,-0.481212,1.1668149,0.83390534,1.0459764,0.9306606,-0.2064376,1.4802107,2.194787,2.1997948,1.4666415,0.4752782,-0.3851213,0.39982617,1.1049168,1.5582439,1.1441491,0.54350895,-0.12377602,-1.1162694,-0.81837773,-0.61954826,-0.5081537,0.05654518,0.27608633,0.058893226,0.66271317,-0.5643443,-1.4979934,-1.6173791,-0.93975836,0.25631756,0.59934616,2.843477,0.89051753,-1.0120623,-2.0241823,-1.226678,1.0279679,-0.8218907,-1.8627697,0.47883475,0.79846036,0.11416366,0.26554617,-0.8187242,-0.67032385,0.022432745,0.12340537,-1.5547149,-1.367483,-0.486928,-0.25378674,-0.72260815,-1.2786175,-0.97312313,-0.6477252,-0.27312654,-0.26127976,-0.72882843,-1.0788345,-0.1731066,0.1412026,0.23371555,-0.039848015,-0.42800424,-0.6013456,-0.9716911,0.08407044,0.12264694,0.07915063,0.068612374,-0.26937574,-0.47231367,-0.85249305,-0.49820036,-0.5340724,-0.51410854,-0.25313598,-0.106181026,-0.15548608,-0.44849333,-0.8423466,-1.1075292,-1.0955443,-0.6117004,0.34099028,-0.1569541,-0.6340152,-1.1396832,-1.7039739,0.07533748,-0.08189709,0.19304663,0.018776672,0.18971433,-0.43763462,-0.37614876,-1.1094376,-0.8876297,-0.24171622,-0.07213691,-1.5315721,-1.3601786,-0.835767,-0.3647269,-0.41614193,-0.7048859,-0.9596979,-0.4280433,-0.2780079,-0.4078675,-0.8567372,-0.71017444,-1.2269157,-0.037848685,0.6420395,0.3552197,-0.4960394,-0.7127392,-0.45230854,0.6883611,1.0982074,0.57119787,-0.24441506,-0.32687584,0.4152601,0.6467641,0.27235496,-0.076504685,-0.6694386,-0.5709675,-0.19175994,-0.07078535,-0.5970133,-0.19487779,-0.33990598,-0.39548707,-0.047516342,0.38796777,0.33274388,-0.101605676,-0.29784444,0.15943407,0.24134225,0.01069601,-0.21617793,0.49520165,0.27437946,-0.20740907,-0.50444955,-0.36690298,1.5124941,0.6566162,-0.5695327,-0.65782535,-0.154533,2.6197684,2.4616811,1.2959715,1.7798659,2.0666533,1.6456947,-0.56362814,-0.5398235,-0.289549,0.18546073,0.73811233,-0.106018014,0.6849294,0.7390298,0.3297052,0.23857288,0.36811852,-0.07911716,-0.2801976,1.1529573,0.44330063,2.2662222"""

output_em="""
1.0551954e-06,0.43045548,1.4624727e-06,1.0827124e-08,1.4222807e-09,0.02086455,6.233141e-10,1.2245289e-12,1.6039878e-10,4.204569e-07,1.6519078e-10,2.3467856e-11,8.918103e-12,6.78244e-08,3.935764e-07,0.5199692,1.2221555e-09,0.0025662058,3.1933683e-10,0.026141122"""

input_xray="""
2.1942575,0.52937645,-1.081864,-1.0748761,0.9040855,1.7355636,0.4674884,-1.3136659,-1.5414683,-0.6928854,0.47842777,-0.043354426,-0.40737396,2.7822642,0.26603904,-0.37640825,-0.43451056,-0.6202667,-0.0043380107,-1.5256929,-0.7986468,-1.2786487,-1.5538476,-1.4740794,-1.2767115,-1.189284,-1.400935,-1.066136,-0.36615312,-1.3615373,-1.4397273,-0.7626196,-0.30828843,2.0490258,-0.846023,-1.7642158,-0.401722,-1.2110116,-0.10815647,-0.9900087,0.5684164,-1.8321794,-1.9446927,-1.9443603,-0.7885904,0.23058867,-1.076702,-1.4191272,-1.2086936,2.0415828,0.060899135,-1.0456253,-0.9088898,0.08076783,0.9622302,-0.65011543,-0.82102424,-0.71820945,-0.4241632,0.13254446,1.0257747,-0.26908413,-0.8361191,-0.3200535,0.082442224,-1.018161,-0.22203428,-0.9980255,-1.186957,-0.45880395,0.083343305,0.057712723,1.0711724,-0.39030492,-1.3429877,-0.5645192,-0.0253975,-0.10595601,2.2662868,-0.15048686,-0.129302,1.5770833,-0.060515128,-0.93629265,-0.71385574,-0.5292911,0.07171272,-0.6656201,0.2648466,0.9151036,-1.2651739,-0.5638953,0.21791226,0.71897334,0.9408783,0.08272985,-0.6496161,-0.53518903,0.35616326,3.1835113,1.8144917,-0.4231017,-0.9183731,-0.39559674,0.6184881,3.3408267,2.6184323,0.6691932,0.2129578,-0.5080888,0.097920194,1.0595902,1.8606575,1.028438,0.82441026,0.2451042,0.732248,-0.7026622,-0.65338486,-0.7985776,-0.20756504,0.16216528,0.40238884,1.2192866,1.0750809,-0.5707935,-1.248426,-0.791049,-0.17870936,-0.43752784,4.363548,3.0712547,0.6669736,-0.53081733,2.830938,0.8517582,-0.56536555,-0.5571253,-1.150689,-0.38207093,0.08277288,-1.5332841,-0.69779813,-0.12526883,-0.8804289,-0.5215831,0.02933978,0.74250585,0.99104893,-0.109344296,-0.59898186,-0.9977784,1.104299,1.5246937,1.6123486,1.129486,-0.86699194,-0.19016097,-0.04413379,0.5426886,1.6928672,1.7783824,0.42406112,0.8281698,-0.7091212,-0.88489145,-1.1037629,-0.21046181,0.64070976,0.5303163,0.68835324,1.5240108,0.891341,-0.23059048,-0.53700805,-0.8497232,-0.43003383,-0.3103036,3.5375373,2.8407757,2.4253936,1.0937947,-0.85214823,1.148648,2.2123911,1.2264115,-0.20163144,0.10368292,-0.10716472,0.05858874,-0.86332536,-0.8289826,-0.6211861,-0.8889359,-1.0025527,-1.0079216,-0.80021065,-0.93996793,-0.052355558,-0.8530387,-0.34673482,-0.2715981,0.015444712,-0.23406316,-0.26459047,-1.0579004,-0.63807744,-0.10639029,1.1829724,1.2226993,0.1620937,0.08601169,-1.4197528,-1.3331912,-0.61691713,0.93203056,1.2714493,0.77362907,1.2046968,-1.4253598,-1.7721952,-1.5324062,-0.30463603,0.16630055,0.28163633,0.87171257,-0.5575599,-1.1910441,-1.2454479,-0.8953918,-1.0316433,-0.795933,0.35045013,-0.48251423,-0.15496409,-0.36221197,0.348512,0.1227249,0.3952529,0.7206121,0.5219623,-0.048127417,0.51832986,-0.5802352,-0.4755818,-0.32925525,-0.7256249,-0.4130007,-1.33436,-1.027688,-0.057869922,-0.023229657,0.35477734,-0.56127745,-1.333919,-0.67185676,0.5643201,0.79507816,0.853372,-0.08494951,-0.6688172,-0.98620605,-0.3885133,0.5769416,0.79580605,1.1729722,0.34812042,-0.33997676,-1.1562874,-0.42815244,0.060444564,1.1151845,0.47878045,-0.78593284,-0.7912496,-0.094627,0.2598724,0.903729,0.3441146,0.37647286,0.8798864,0.78296924,0.6244839,0.3991105,0.27588606,0.08486053,0.1337329,-0.11789984,-0.09011547,0.5511347,0.37435585,0.22735311,0.08188228,0.035412405,0.1655578,-0.011522835,0.40467787,-0.26199138,-0.15851824,0.026246678,0.103138894,0.74468386,0.95235246,0.26887876,-0.5862246,1.5990161,1.108012,-0.010009038,0.30963427,0.504617,0.22236821,0.34292325,0.3095476,-0.051055346,-0.21903951,0.32758477,0.5363083,0.12710758,-0.2790137,-0.21343572,-0.045869797,0.027584942,0.4022605,0.045373715"""

output_xray="""
1.5994388e-10,0.00083333603,1.248885e-09,1.3974322e-06,1.12567206e-10,0.96025276,1.6565794e-19,7.441627e-19,9.167633e-15,6.260032e-10,7.940784e-11,1.4934677e-14,1.1789899e-14,5.7310075e-08,1.3321625e-07,0.031721976,2.9224486e-06,0.0071872906,1.2573872e-12,6.690157e-10"""



ROOT = os.path.dirname(os.path.realpath(__file__))


import torch

from torch.autograd import Variable
import torch.nn.functional as F
from torch import optim
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset

from fmslib import nn_models

import numpy as np

from fmslib import descriptors


from torchsummary import summary

def main(em=1):

    if em:
        print(" ==> Parsed EM model")
        fname=os.path.join(ROOT, '..', 'data', 'nn_gpu_1e4_bs200_naive.dat')
        tscript_fname=os.path.join(ROOT, 'em_model.pt')
    else:
        print(" ==> Parsed XRAY model")
        fname=os.path.join(ROOT, '..', 'data', 'xtst_5k.dat')
        tscript_fname=os.path.join(ROOT, 'xray_model.pt')

    device = torch.device("cpu")

    model = nn_models.naive().to(device)

    model.load_state_dict(torch.load(fname, map_location=device))
    model.eval()

    if em:
        input_np = np.array(input_em.split(','), dtype=np.float32)
        output_np = np.array(output_em.split(','), dtype=np.float32)
    else:
        input_np = np.array(input_xray.split(','), dtype=np.float32)
        output_np = np.array(output_xray.split(','), dtype=np.float32)


    input_var = Variable(torch.FloatTensor(input_np))
    input_var = torch.tensor(np.array([input_np])).float()

    print(input_var.shape)
    err = np.linalg.norm(output_np-model(torch.tensor([input_np]).float()).exp().detach().numpy().flatten())
    assert err<1e-10


    #k2c_input_np=np.array(k2c_input_str.replace('f','').split(","), dtype=np.float32)
    #k2c_output_np=np.array(k2c_output_str.replace('f','').split(","), dtype=np.float32)

    #_m = model(torch.tensor([k2c_input_np]).float()).exp().detach().numpy().flatten()

    #print(_m-np.exp(k2c_output_np)/np.sum(np.exp(k2c_output_np)))

    torch.onnx.export(model, input_var, "%s_model.onnx" % ('em' if em else 'xray') )
    summary(model, (1,324))


if __name__=="__main__":
    main(em=1)
    print('-')
    main(em=0)
