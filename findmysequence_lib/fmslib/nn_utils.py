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



import sys, os, re

import h5py
import numpy as np
from fmslib import xyz_utils
from fmslib import descriptors

class nn_utils:

    def __init__(self):

        self._xyz_utils = xyz_utils.xyz_utils()
        pass

    # -------------------------------------------------------------------------

    def write_training_set(self, modelin,
                                 mapin      =   None,
                                 mtzin      =   None,
                                 labin      =   None,
                                 h5py_fname =   None):

        _descr = descriptors.descriptors(mapin=mapin, mtzin=mtzin, labin=labin)
        ph, symm = _descr._xyz_utils.read_ph(modelin)
        desc_data = {}

        if h5py_fname is None:
            h5py_fname = 'nn_trainset_%s.hdf5'%( os.path.splitext(os.path.basename(modelin))[0] )

        for ires,res in enumerate(ph.atom_groups()):
            res_type, cloud = _descr.describe(res)
            if cloud is None: continue

            _a = desc_data.setdefault(res_type, [])
            _a.append(cloud)


        with h5py.File(h5py_fname, 'w') as data_file:
            for key in desc_data.keys():
                print( key, len(desc_data[key]) )
                data_file.create_dataset(key, data=np.array(desc_data[key]))


        return os.path.abspath(h5py_fname)


    # -------------------------------------------------------------------------

    def read_training_set(self, h5py_fname):
        """
            data - array of 
        """

        data, labels = [], []
        with h5py.File(h5py_fname, 'r') as data_file:
            for _k, _v in data_file.items():
                # there should be only std names, but...
                if not _k in self._xyz_utils.standard_aa_names: continue

                data.extend(_v)
                labels.extend( np.ones(len(_v))*self._xyz_utils.standard_aa_names.index(_k) )

        assert len(labels) == len(data)

        return data, labels#np.array(data), np.array(labels)


    # --------------------------

    def read_multiple_sets(self, fnames):


        data, labels = [], []

        for fn in fnames:
            d,l = self.read_training_set(fn)
            print(fn, len(d))
            data.extend(d)
            labels.extend(l)

        data = np.array(data)
        labels=np.array(labels)

        return data, labels

    # --------------------------

    def select_trainset(self, data, labels, valid_fraction=0.1):

        n = len(data)
        # select test_set
        idx = np.random.choice(n, int(valid_fraction*n), replace=False)

        train_data, train_labels = np.delete(data, idx, axis=0), np.delete(labels, idx)
        test_data, test_labels = data[idx,:], labels[idx]

        assert len(labels) == (len(test_data) + len(train_data))

        return train_data, train_labels, test_data, test_labels


