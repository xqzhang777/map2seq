#! /use/bin/env libtbx.python

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

import os

from cctbx import array_family
from fmslib import xyz_utils

# scitbx
import scitbx
from scitbx import matrix
from scitbx.math import superpose
import numpy as np


ROOT = os.path.dirname(os.path.realpath(__file__))

class descriptors:



    def __init__(self, mapin=None, mtzin=None, labin=None, badd=None):

        self._xyz_utils = xyz_utils.xyz_utils()
        self.bb_ph, self.cloud_ph = self._xyz_utils.read_cloud(ifname = \
                                            os.path.join(ROOT, '..', 'data', 'cloud.pdb'))

        self.template_bb_xyz = array_family.flex.vec3_double(self._xyz_utils.res_bb(self.bb_ph.only_residue_group()))

        # map has a preference
        self._xyz_utils.init_maps(mapin=mapin, mtzin=mtzin, labin=labin, badd=badd)


    # -------------------------------------------------------------------------



    def describe(self, res, strict=False):


        try:
            _xyz = array_family.flex.vec3_double(self._xyz_utils.res_bb(res))
        except:
            return None, None

        res_type = res.resname.strip().upper()

        if strict and not res_type in self._xyz_utils.standard_aa_names:
            return None, None

        superposition = superpose.least_squares_fit(_xyz, self.template_bb_xyz, \
                                                    method=["kearsley", "kabsch"][0])

        rtmx = matrix.rt((superposition.r, superposition.t))

        cloud = self.cloud_ph.only_residue_group().detached_copy()

        cloud.atoms().set_xyz( rtmx * cloud.atoms().extract_xyz() )
        #return cloud
        _cloud = []
        for atm in sorted(cloud.atoms(), key=lambda _a: _a.i_seq):
            map_val = self._xyz_utils.map_data.eight_point_interpolation(
                          self._xyz_utils.uc.fractionalize(atm.xyz))
            _cloud.append(map_val)

        _cloud = np.array(_cloud)
        std,mean = _cloud.std(), _cloud.mean()
        if std>0: _cloud = (_cloud-mean)/std

        return res_type, _cloud


def tst2():
    import iotbx
    tmp_frag = iotbx.pdb.hierarchy.root()
    tmp_frag.append_model(iotbx.pdb.hierarchy.model(id="0"))
    tmp_frag.models()[0].append_chain(iotbx.pdb.hierarchy.chain(id="0"))

    _descr = descriptors()
    ph, symm = _descr._xyz_utils.read_ph('helix.pdb')
    for ires,res in enumerate(ph.atom_groups()):
        cloud = _descr.describe(res)
        if ires%5==0:
            tmp_frag.only_chain().append_residue_group( cloud.detached_copy() )

    with open("helical_cloud.pdb", 'w') as of:
        of.write(tmp_frag.as_pdb_string())

def tst(modelin='../examples/5me2.pdb', mapin='../examples/emd_3488.map'):
    _descr = descriptors(mapin)
    ph, symm = _descr._xyz_utils.read_ph(modelin)
    for ires,res in enumerate(ph.atom_groups()):
        _t, _c = _descr.describe(res)
        if _c is None: continue
        print("%5i %5s %3.2f" % (ires, _t, _c.std()))
        if ires>10: break


if __name__=="__main__":
    tst()
