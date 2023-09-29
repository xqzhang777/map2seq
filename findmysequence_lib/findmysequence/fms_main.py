#! /usr/bin/env python

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
__date__ = "08 Apr 2021"


import os,sys,re
from shutil import which

#from optparse import OptionParser, OptionGroup, SUPPRESS_HELP

ROOT = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(2, '%s/../fmslib' % ROOT)
sys.path.insert(2, '%s/../' % ROOT)

import numpy as np

import fmslib
from fmslib import xyz_utils
from fmslib import sequence_utils

HMMER_AVAILABLE = which('hmmsearch')

try:
    import findmysequence.version
    version = findmysequence.version.__version__
except:
    version="dev"

fms_header="""
# findMySequence(v%s): tool for finding and assigning protein model sequences in EM and MX
"""%(version)

class Fms_option:
    def __init__(self,mapin,modelin,seqin, modelout, slide, db, tmpdir,outdir,rev,flip,tophits=3):
        self.mapin=mapin
        self.modelin=modelin
        self.seqin=seqin
        self.modelout=modelout
        self.db=db
        self.tmpdir=tmpdir
        self.outdir=outdir
        self.tophits=tophits
        #self.slide=True if HMMER_AVAILABLE is None else False
        self.slide=slide
        self.selstr="all"
        self.rev=int(rev)
        self.flip=int(flip)

#def parse_args():
#    """setup program options parsing"""
#    parser = OptionParser()


#    required_opts = OptionGroup(parser, "Required parameters (model and a map or an mtz file with labels)")
#    parser.add_option_group(required_opts)

#    required_opts.add_option("--modelin", action="store", \
#                            dest="modelin", type="string", metavar="FILENAME", \
#                  help="PDB/mmCIF model", default=None)

#    required_opts.add_option("--mapin", action="store", dest="mapin", type="string", metavar="FILENAME", \
#                  help="Map in CCP4 format", default=None)

#    required_opts.add_option("--mtzin", action="store", dest="mtzin", type="string", metavar="FILENAME", \
#                  help="input sfs in MTZ format", default=None)

#    required_opts.add_option("--labin", action="store", dest="labin", type="string", metavar="F,PHI", \
#                  help="MTZ file column names", default=None)
    
#    required_opts.add_option("--tmpdir", action="store", dest="tempdir", type="string", default=None)


#    extra_opts = OptionGroup(parser, "Extra options")
#    parser.add_option_group(extra_opts)

#    extra_opts.add_option("--select", action="store", dest="selstr", type="string", metavar="STRING", \
#                  help="fragments selection string"\
#                       " eg. \"resi 10:50 and chain A\" [default: %default]", default="all")

#    extra_opts.add_option("--tophits", dest="tophits_sto", type="int", default=3, metavar="VALUE", \
#                  help="Number of top sequence hits to print [default: %default]")

#    default_ccp4_db=None
#    for fn in ['pdb100.txt', 'pdbALL.txt']:
#        _default_ccp4_db = os.path.expandvars("$CCP4/share/mrbump/data/%s"%fn)
#        if os.path.exists(_default_ccp4_db):
#            default_ccp4_db = _default_ccp4_db
#            break

#    ccp4 = os.environ.get('CCP4', None)

#    # DB-serach related options
#    search_opts = OptionGroup(parser, "Database search options (default mode)")
#    parser.add_option_group(search_opts)

#    search_opts.add_option("--slide", action="store_true", dest="slide", default=True if HMMER_AVAILABLE is None else False, \
#            help="Use naive sliding-window algorithm for db search (much slower than default hmmer) [default: %default]")

#    search_opts.add_option("--db", action="store", \
#                            dest="db", type="string", metavar="FILENAME", \
#                            help="Sequence database (may be gzipped) [default: %default]", default=default_ccp4_db)


#    # sequende assignment specific options
#    dock_opts = OptionGroup(parser, "Sequence assignement options (both arguments below are required)")
#    parser.add_option_group(dock_opts)


#    dock_opts.add_option("--seqin", action="store", dest="seqin_fname", type="string", metavar="FILENAME", \
#                  help="targt sequence for sequene assignment", default=None)

#    dock_opts.add_option("--modelout", action="store", dest="modelout", type="string", metavar="FILENAME", \
#                  help="output filename for the sequence assignment", default=None)


#    # developer opts
#    parser.add_option("--debug", action="store_true", dest="debug", default=False, \
#                  help=SUPPRESS_HELP)

#    parser.add_option("--refseq", action="store", dest="refseq_fname", type="string", metavar="FILENAME", \
#                  help=SUPPRESS_HELP if 1 else """reference sequence in fasta format (for testing, will print 
#                                                    sequence identity to all hits)""", default=None)

#    parser.add_option("--test", action="store_true", dest="test", default=False, \
#                  help="a simple test")

#    required_opts.add_option("--badd", action="store", dest="badd", type="float", metavar="FLOAT", \
#                  help=SUPPRESS_HELP, default=0)


#    (options, _args)  = parser.parse_args()
#    return (parser, options)


# -----------------------------------------------------------------------------


def guess(mapin         =   None,
          mtzin         =   None,
          labin         =   None,
          modelin       =   None,
          db            =   None,
          selstr        =   None,
          refseq_fname  =   None,
          modelout      =   None,
          slide         =   False,
          tophits_sto   =   3,
          badd          =   0,
          debug         =   False,
          tempdir       =   None,
          outdir        =   None,
          rev           =   0,
          flip          =   0):

    m2so = sequence_utils.model2sequence( mapin=mapin, mtzin=mtzin, labin=labin, badd=badd)
    ph, symm = m2so._xyz_utils.read_ph( modelin )
    msa_string = m2so.model2msa( ph, selstr=selstr, verbose=False )
    if rev == 1:
        print("reversing msa")
        msa_string=reverse_msa(msa_string)


    print(" ==> Querying database using %s: %s" % ('sliding window' if slide else 'hmmer', db))
    if slide:
        res = m2so.query_slide(db_fname=db, refseq_fname=refseq_fname, verbose=debug, tophits_sto=tophits_sto)
    else:
        res = m2so.query_msa( msa_string, db_fname=db, refseq_fname=refseq_fname, verbose=debug, tophits_sto=tophits_sto)
    stdout = sys.stdout
    sys.stdout = open(f'{outdir}/hmmer_output.txt', 'w')
    if res is None:
        print(" ==> No matches found...")
    else:
        print(" ==> Best matches")
        for v in res:
            print( v[0], v[1] )
        sys.stdout = stdout
        print("Finished")
    save_residue_score_dict(m2so.residue_scores_dicts,outdir)
    return res


def save_residue_score_dict(score_dict,outdir):
    import pandas as pd
    import pickle  

    df=pd.DataFrame.from_dict(score_dict,orient="index")
    with open(f'{outdir}/score_dict.pkl','wb') as o:
        pickle.dump(df,o,pickle.HIGHEST_PROTOCOL)
# -----------------------------------------------------------------------------

def assign_sequence(mapin         =   None,
                    mtzin         =   None,
                    labin         =   None,
                    modelin       =   None,
                    db            =   None,
                    selstr        =   None,
                    seqin_fname   =   None,
                    modelout      =   None,
                    slide         =   False,
                    debug         =   False,
                    outdir        =   None,
                    rev           =   0):

    m2so = sequence_utils.model2sequence( mapin=mapin, mtzin=mtzin, labin=labin )
    ph, symm = m2so._xyz_utils.read_ph( modelin )
    msa_string = m2so.model2msa( ph, selstr=selstr, verbose=False )
    if rev == 1:
        print("reversing msa")
        msa_string=reverse_msa(msa_string)
    
    stdout = sys.stdout
    sys.stdout = open(f'{outdir}/seq_align_output.txt', 'w')
    print( " ==> Aligning chain fragments to the input sequence" )
    with open(seqin_fname, 'r') as ifile:
        m2so.align_frags(target_sequence=ifile.read(), modelout=modelout, verbose=debug)
    sys.stdout = stdout
    print("Finished")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def fms_run(mapin=None, modelin=None, seqin=None, modelout=None, slide=False, db=None, tmpdir=None, outdir=None, tophits=None, rev=0, flip=0):

#    (parser, options) = parse_args()

#    print( " ==> Command line: %s" % (" ".join(sys.argv)) )

#    if options.test:
        
#        print("\n")
#        print( " *** Testing EM search\n")
#        guess( mapin    = os.path.join(xyz_utils.ROOT, '..', 'examples', 'emd_3488.map'),
#               modelin  = os.path.join(xyz_utils.ROOT, '..', 'examples', '5me2.pdb'),
#               selstr   = "chain A and resi 20:40",
#               slide    = True,
#               debug    = options.debug,
#               db       = os.path.join(xyz_utils.ROOT, '..', 'examples', 'example_db.fa.gz'))

#        print("\n")
#        print( " *** Testing XRAY search")
#        guess( mtzin    = os.path.join(xyz_utils.ROOT, '..', 'examples', '1cbs_final.mtz'),
#               labin    = "FWT,PHWT",
#               modelin  = os.path.join(xyz_utils.ROOT, '..', 'examples', '1cbs_final.pdb'),
#               selstr   = "chain A and resi 20:40",
#               slide    = True if HMMER_AVAILABLE is None else False,
#               debug    = options.debug,
#               db       = os.path.join(xyz_utils.ROOT, '..', 'examples', 'example_db.fa.gz'))


#        print("\n")
#        print( " *** Testing sequence assignment in EM")
#        assign_sequence( mapin        =   os.path.join(xyz_utils.ROOT, '..', 'examples', 'emd_3488.map'),
#                         modelin      =   os.path.join(xyz_utils.ROOT, '..', 'examples', '5me2.pdb'),
#                         selstr       =   "chain A and resi :20",
#                         seqin_fname  =   os.path.join(xyz_utils.ROOT, '..', 'examples', '5me2.fa'),
#                         modelout     =   None)

#        exit(0)

#    valid_params_no = len(sys.argv[1:]) - len(parser.largs)

    ## no recognized params on input, print help message and exit...
#    if not valid_params_no:
#        parser.print_help()
#        print
#        return 1
    sld123123=slide
    options=Fms_option(mapin=mapin, modelin=modelin, seqin=seqin, modelout=modelout, db=db, tmpdir=tmpdir, outdir=outdir, tophits=tophits,rev=rev,flip=flip, slide=sld123123)

#    if options.modelin is None:

#        print("ERROR: Input model missing")
#        print()
#        return 1

#    if options.mapin is None and \
#        (options.mtzin is None and options.labin is None):

#        print("ERROR: Input map or mtz file with labels missing")
#        print()
#        return 1


#    if [options.seqin_fname,options.modelout].count(None)==1:
#        print("ERROR: Both seqin and modelout are required for the sequence asignment")
#        print()
#        return 1

#    if options.seqin_fname and options.modelout:
#        assign_sequence( mapin        =   options.mapin,
#                         mtzin        =   options.mtzin,
#                         labin        =   options.labin,
#                         modelin      =   options.modelin,
#                         db           =   options.db,
#                         selstr       =   options.selstr,
#                         seqin_fname  =   options.seqin_fname,
#                         modelout     =   options.modelout,
#                         slide        =   options.slide,
#                         debug        =   options.debug)
#        return 0

    if options.seqin and options.modelout:
        assign_sequence( mapin        =   options.mapin,
                         mtzin        =   None,
                         labin        =   None,
                         modelin      =   options.modelin,
                         db           =   options.db,
                         selstr       =   options.selstr,
                         seqin_fname  =   options.seqin,
                         modelout     =   options.modelout,
                         slide        =   options.slide,
                         debug        =   False,
                         outdir       =   options.outdir,
                         rev          =   options.rev)
        return 0


    if options.db is None or not os.path.exists(options.db):
        print("ERROR: Sequence database not available")
        print()
        return 1


#    guess( mapin        =   options.mapin,
#           mtzin        =   options.mtzin,
#           labin        =   options.labin,
#           modelin      =   options.modelin,
#           db           =   options.db,
#           selstr       =   options.selstr,
#           refseq_fname =   options.refseq_fname,
#           modelout     =   options.modelout,
#           slide        =   options.slide,
#           tophits_sto  =   options.tophits_sto,
#           badd         =   options.badd,
#           debug        =   options.debug,
#           tempdir       =  options.tempdir)

    hmm_res=guess( mapin        =   options.mapin,
           mtzin        =   None,
           labin        =   None,
           modelin      =   options.modelin,
           db           =   options.db,
           selstr       =   options.selstr,
           refseq_fname =   None,
           modelout     =   None,
           slide        =   options.slide,
           tophits_sto  =   options.tophits,
           badd         =   0,
           debug        =   False,
           tempdir      =   options.tmpdir,
           outdir       =   options.outdir,
           rev          =   options.rev,
           flip         =   options.flip)
    return hmm_res

def reverse_msa(msa_string):
    lines=msa_string.split()
    new_msa=""
    for line in lines:
        if line[0] == ">":
            new_msa+=line+'\n'
        else:
            new_msa+=line[::-1]+'\n'
    return new_msa


if __name__=="__main__":
    fms_run()
