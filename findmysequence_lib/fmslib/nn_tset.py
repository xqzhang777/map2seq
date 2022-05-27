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

ROOT = os.path.dirname(os.path.realpath(__file__))

import numpy as np

from nn_utils import nn_utils
from xyz_utils import xyz_utils


from optparse import OptionParser, SUPPRESS_HELP

# -----------------------------------------------------------------------------

def parse_args():
    """setup program options parsing"""
    parser = OptionParser()


    parser.add_option("--modelin", action="store", \
                            dest="modelin", type="string", metavar="FILENAME", \
                  help="Input PDB/mmCIF model", default=None)

    parser.add_option("--mapin", action="store", dest="mapin", type="string", metavar="FILENAME", \
                  help="input map in CCP4 format", default=None)

    parser.add_option("--mtzin", action="store", dest="mtzin", type="string", metavar="FILENAME", \
                  help="input sfs in MTZ format", default=None)

    parser.add_option("--labin", action="store", dest="labin", type="string", metavar="F,PHI,[FOM]", \
                  help="MTZ file column names", default=None)

    parser.add_option("-o", action="store", dest="ofname", type="string", metavar="FILENAME", \
                  help="output hdf5 filename", default=None)

    parser.add_option("--test", action="store_true", dest="test", default=False, \
                  help="run a basic test")

    (options, _args)  = parser.parse_args()
    return (parser, options)



# -----------------------------------------------------------------------------

def basic_test():
    _nn_utils = nn_utils()
    fn = _nn_utils.write_training_set(mapin   = os.path.join(ROOT, '..', 'examples', 'emd_3488.map'),
                                      modelin = os.path.join(ROOT, '..', 'examples', '5me2.pdb'),
                                      h5py_fname='tst.hdf5')

    fn = _nn_utils.write_training_set(mtzin   = os.path.join(ROOT, '..', 'examples', '1cbs_final.mtz'),
                                      labin   = "FWT,PHWT",
                                      modelin = os.path.join(ROOT, '..', 'examples', '1cbs_final.pdb'),
                                      h5py_fname= 'tst.hdf5')



def create_trainig_set(mapin=None, modelin=None, mtzin=None, labin=None, ofname=None):

    _nn_utils = nn_utils()
    fn = _nn_utils.write_training_set(mapin      = mapin,
                                      mtzin      = mtzin,
                                      labin      = labin,
                                      modelin    = modelin,
                                      h5py_fname = ofname)
    print("Wrote %s" % fn)



# -----------------------------------------------------------------------------

def main():

    (parser, options) = parse_args()

    print( " ==> Command line: %s" % (" ".join(sys.argv)) )

    if options.test:
        basic_test()
        exit(0)

    if options.modelin is None or \
        (options.mapin is None and None in [options.mtzin, options.labin]):
        parser.print_help()
        exit(1)

    create_trainig_set(mapin   = options.mapin,
                       modelin = options.modelin,
                       mtzin   = options.mtzin,
                       labin   = options.labin,
                       ofname  = options.ofname)




if __name__=="__main__":
    main()
