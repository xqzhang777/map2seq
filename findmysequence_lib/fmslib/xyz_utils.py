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
__date__ = "16 Feb 2020"


import sys, os, re

# it seems that real_space_refinement_simple
# supports greedy multiprocessing, that for simple,
# and fast single-sc fitting makes no good
os.environ["OMP_NUM_THREADS"] = "1"



ROOT = os.path.dirname(os.path.realpath(__file__))



import random
import subprocess
import string

import numpy as np

# cctbx imports

from cctbx  import xray
from cctbx  import maptbx
from cctbx  import miller

from cctbx.array_family import flex
from cctbx import crystal
from cctbx import sgtbx, uctbx

import iotbx
import iotbx.map_tools
from iotbx import crystal_symmetry_from_any
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.pdb import amino_acid_codes as aac
import iotbx.pdb
import iotbx.ccp4_map

from mmtbx.maps import utils

# mmtbx
from mmtbx.monomer_library import idealized_aa
from mmtbx.monomer_library import server
from mmtbx.geometry_restraints import reference
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager
from cctbx.maptbx import real_space_refinement_simple

from libtbx.math_utils import ifloor, iceil

# scitbx
import scitbx
from scitbx import matrix
from scitbx import fftpack
from scitbx.math import superpose


import tempfile

standard_aa_names = "CYS,GLN,ILE,SER,VAL,LYS,TRP,PRO,GLY,THR,ALA,PHE,HIS,ASN,LEU,ARG,ASP,GLU,TYR,MET".split(",")

class xyz_utils:


    def __init__(self, logfile=None):

        self.logfile = logfile


        self.standard_aa_names = standard_aa_names
        self.standard_na_names = "A,G,U,C".split(",")

        self.ncaco_dict = {'N':0, 'CA':1, 'C':2}

        self.alphanumbers = "ABCDEFGHIJKLMNOPabcdefghijklmnopqrstuvwxyz0123456789RSTUVWXYZ"

        self.workdir = '.'

        self.ogt = aac.one_letter_given_three_letter
        self.tgo = aac.three_letter_given_one_letter

        self.ideal_dict = idealized_aa.residue_dict()

        if not "CLIBD_MON" in os.environ: os.environ["CLIBD_MON"] = os.path.join(ROOT, "..", "data", "mon_lib")
        self.mon_lib_srv = server.server()

        # populate mon_lib_srv with local cif data (files provided by CCP4 doesn't work)
        # here, standard residues are enough
        # source: https://sourceforge.net/projects/geostd/

        for _rn in self.standard_aa_names:
            self.mon_lib_srv.process_cif( file_name=os.path.join(ROOT, "..", "data", "geostd", _rn[0].lower(), "data_%s.cif"%_rn.upper()) )


        self.dmin = 3.5

    # -------------------------------------------------------------------------

    def _logger(self, msg):

        if self.logfile:
            self.logfile.write("%s\n" % msg)
        else:
            print( msg )

    # -------------------------------------------------------------------------

    def res_bb(self, res):

        ncaco_xyz = [None]*len(self.ncaco_dict.keys())

        for atom in res.atoms():
            name = atom.name.strip()
            try:
                ncaco_xyz[ self.ncaco_dict[name] ] = atom.xyz
            except:
                pass

        return ncaco_xyz


    # -------------------------------------------------------------------------

    def parse_pdbstring(self, pdb_string):


        # there may be issues with repeated BREAK lines, that we do not use here anyway
        pdb_string_lines = []
        for _line in pdb_string.splitlines():
            if re.match(r'^BREAK$', _line):
                continue
            pdb_string_lines.append(_line)


        # arghhhh, why do you guys keep changing the interface?
        inp = iotbx.pdb.input(source_info=None, lines=pdb_string_lines)
        try:
            return inp.construct_hierarchy(sort_atoms=False), inp.crystal_symmetry()
        except:
            return inp.construct_hierarchy(), inp.crystal_symmetry()


    # -------------------------------------------------------------------------

    def read_ph(self, ifname, verbose=True):

        if verbose: self._logger( " ==> Parsing a PDB/mmCIF file: %s" % ifname )


        with open(ifname, 'r') as ifile:
            ph, symm = self.parse_pdbstring(ifile.read())


        ph.remove_alt_confs(True)

        return ph, symm


    # -------------------------------------------------------------------------


    def find_chain_fragments(self, ph, verbose=True, fake=0):


        acepted_aa_names = self.standard_aa_names

        chain_fragments = []

        # find fragments of chains with iseq gaps in between
        # TODO: something more sophisticated maye be needed
        for chain in ph.chains():

            if verbose: self._logger( " ==> Processing chunks in chain %s" % chain.id )

            res_chunk = []
            for rg in chain.residue_groups():

                # check if CA is there, if not and natoms is larger that 1 it must be aa
                if not "CA" in [_.name.strip() for _ in rg.atoms()] or rg.atoms().size()==1: continue


                resname = rg.only_atom_group().resname
                resseq = rg.resseq_as_int()

                # ignore dums, in case input from arpwarp is used
                if resname == "DUM":
                    continue

                if res_chunk and res_chunk[-1].resseq_as_int()<(resseq-1):
                    chain_fragments.append(res_chunk)
                    res_chunk = [rg]
                    continue

                res_chunk.append(rg)

            if res_chunk:
                chain_fragments.append(res_chunk)


        return chain_fragments

    # -------------------------------------------------------------------------


    def find_chain_fragments_na(self, ph, verbose=True, fake=0):


        acepted_resi_names = self.standard_na_names

        chain_fragments = []

        # find fragments of chains with iseq gaps in between
        # TODO: something more sophisticated maye be needed
        for chain in ph.chains():

            if verbose: self._logger( " ==> Processing chunks in chain %s" % chain.id )

            res_chunk = []
            for rg in sorted(chain.residue_groups(), key=lambda _: _.resseq_as_int()):

                # check if CA is there, if not and natoms is larger that 1 it must be aa
                if not "P" in [_.name.strip() for _ in rg.atoms()] or rg.atoms().size()==1: continue


                # get rid of alternative conformations
                if len(rg.atom_groups())>1:
                    for _ag in rg.atom_groups()[1:]: rg.remove_atom_group(_ag)

                resname = rg.only_atom_group().resname
                resseq = rg.resseq_as_int()

                # ignore dums, in case input from arpwarp is used
                if resname == "DUM":
                    continue

                if res_chunk and res_chunk[-1].resseq_as_int()<(resseq-1):
                    chain_fragments.append(res_chunk)
                    res_chunk = [rg]
                    continue

                res_chunk.append(rg)

            if res_chunk:
                chain_fragments.append(res_chunk)


        return chain_fragments

    # -------------------------------------------------------------------------


    # expand input map to cover whole UC
    # cctbx ocasionally has problems with implicit symmetry expansion
    def prepare_input_map(self, mapin, tmpdir=None, verbose=True):
        #if verbose: print " ==> Preparing input map"

        if tmpdir is None: tmpdir = self.workdir

        imap_abspath = os.path.abspath(mapin)
        imap_basename = "".join( os.path.basename(imap_abspath).split('.')[:-1] )

        map_exp = "%(tmpdir)s/%(imap_basename)s_expanded.map" % locals()

        # get input map params
        m = iotbx.ccp4_map.map_reader(file_name=imap_abspath)
        symm_int = m.space_group_number
        epipe = subprocess.Popen("""$CBIN/mapmask mapin %(mapin)s mapout %(map_exp)s << EOF
                                    mode mapin
                                    xyzlim 0 1 0 1 0 1
                                    symm %(symm_int)i
                                    pad 0.0
                                    axis X Y Z
                                    end
                                    EOF""" % locals(), \
                                             shell=True, \
                                             stdout=subprocess.PIPE, \
                                             stderr=subprocess.PIPE)

        retcode = subprocess.Popen.wait(epipe)

        return retcode, map_exp


    # -------------------------------------------------------------------------

    def init_maps(self, mtzin=None, labin=None, mapin=None, debug=False, badd=0, verbose=True):


        assert [mtzin, mapin].count(None)==1

        if mtzin is None:
            self.workdir =  os.path.dirname(os.path.abspath(mapin))
        else:
            self.workdir =  os.path.dirname(os.path.abspath(mtzin))


        if mapin is not None:



            self._logger( " ==> Parsing input CCP4 map from %s" % mapin )

            with tempfile.TemporaryDirectory(prefix="guessmysequence_") as tmpdirname:
                ret, mapin_exp = self.prepare_input_map(mapin, tmpdir=tmpdirname)
                #if mapmask fails (e.g. when CCP4 not available) try to use input maps directly
                if ret: mapin_exp = mapin
                inp_ccp4_map = iotbx.ccp4_map.map_reader(file_name=mapin_exp)


            self.symm = inp_ccp4_map.crystal_symmetry()

            #seems, that this UC may differ from params parsed from the same map e.g. by COOT
            #self.uc = inp_ccp4_map.unit_cell()
            #this is (usually) correct
            self.uc = uctbx.unit_cell(inp_ccp4_map.unit_cell_parameters)

            self.map_data = inp_ccp4_map.data.as_double()

            self.asu_mappings=self.symm.asu_mappings(buffer_thickness=0)






            special_position_settings = crystal.special_position_settings(
                                                    crystal_symmetry=self.symm, \
                                                    min_distance_sym_equiv=0.05)




            n_real = inp_ccp4_map.unit_cell_grid
            grid2frac = maptbx.grid2frac(n_real)
            frac2grid = maptbx.frac2grid(n_real)

            print(" ==> Map is origin shifted: ", not self.map_data.is_0_based())

            if 0:
                s = maptbx.statistics(self.map_data)
                o = self.map_data.origin()
                f = self.map_data.focus()

                self.map_data = maptbx.copy((self.map_data-s.mean())/s.sigma(),
                                            tuple(o), tuple((f[0]-2, f[1]-2, f[2]-2))).as_double()

            if not self.map_data.is_0_based():#not flex.grid(n_real).size_1d() == self.map_data.size():
                # expand the map to the unit cell grid
                # to avoid problems with different ASU definitions
                # and for faster asu-mappings
                map_new = flex.double(flex.grid(n_real), 0)
                o = self.map_data.origin()
                f = self.map_data.focus()
                for i in range(o[0],f[0]):
                    for j in range(o[1],f[1]):
                        for k in range(o[2],f[2]):
                            map_new[i%n_real[0], j%n_real[1], k%n_real[2]] = self.map_data[i, j, k]

                self._logger( "WARNING: problems with map procesing. Is:        %s" % str(self.map_data.accessor().focus()) )
                self._logger( "WARNING: problems with map procesing. Should be: %s" % str(map_new.accessor().focus()) )

                self.map_data = map_new

                if False:
                    iotbx.ccp4_map.write_ccp4_map(
                        file_name="expanded.map",
                        unit_cell=self.uc,
                        space_group=self.symm.space_group_info().group(),
                        gridding_first=self.map_data.accessor().origin(),
                        gridding_last=self.map_data.accessor().last(),
                        map_data=self.map_data,
                        labels=flex.std_string(["iotbx.ccp4_map.tst"]))


            if badd is not None:
                self._logger( " ==> Processing input map (badd %.2fA^2)" % badd )

                crystal_gridding = maptbx.crystal_gridding(
                                               unit_cell             = self.symm.unit_cell(),
                                               space_group_info      = sgtbx.space_group_info(number=1),
                                               pre_determined_n_real = self.map_data.all())


                from iotbx.map_manager import map_manager

                mmgr = map_manager(map_data = self.map_data, unit_cell_grid = self.map_data.all(), unit_cell_crystal_symmetry = self.symm, wrapping=False)

                # max resol 2.0 should be fine, oder?
                map_coeffs = mmgr.map_as_fourier_coefficients(d_min=2.5, d_max=99.99)

                if not map_coeffs:
                    self._logger( "WARNING: map sharpening failed" )
                    return

                #mtz_dataset = map_coeffs.as_mtz_dataset(column_root_label="F")

                d_star_sq = map_coeffs.unit_cell().d_star_sq(map_coeffs.indices())
                dw = flex.exp(-0.25*d_star_sq*badd)

                map_coeffs_sharpen = miller.array(miller_set=map_coeffs, data=map_coeffs.amplitudes().data()*dw)
                map_coeffs_sharpen = map_coeffs_sharpen.phase_transfer(map_coeffs)

                target_map = map_coeffs_sharpen.fft_map(resolution_factor = 1.0/3.0, crystal_gridding = crystal_gridding, symmetry_flags = maptbx.use_space_group_symmetry)
                target_map.apply_sigma_scaling()


                self.map_data = target_map.real_map_unpadded()


                if False:
                    iotbx.ccp4_map.write_ccp4_map(
                    file_name="blurred_generic.map",
                    unit_cell=self.uc,
                    space_group=self.symm.space_group_info().group(),
                    map_data=self.map_data,
                    labels=flex.std_string(["iotbx.ccp4_map.tst"]))


            s = maptbx.statistics(self.map_data)

            self._logger( "     Space group: %s" % self.symm.space_group_info().symbol_and_number() )
            self._logger( "     Unit cell:   %s" % ",".join( [ "%.2f" % x for x in self.uc.parameters()] ) )
            self._logger( "     mean/sigma:  %.2f/%.2f" % (s.mean(),s.sigma()) )
            self._logger( "" )

            return



        self._logger( " ==> Reading input map from MTZ file (%s)" % labin )


        reflection_file = reflection_file_reader.any_reflection_file( file_name = mtzin )
        self.symm = reflection_file.file_content().crystals()[0].crystal_symmetry()

        xray_data_server =  reflection_file_utils.reflection_file_server(
                                                crystal_symmetry = self.symm,
                                                force_symmetry = True,
                                                reflection_files=[reflection_file])



        self.maxres = reflection_file.file_content().max_min_resolution()[1]

        all_labels = []

        for label in utils.get_map_coeff_labels(xray_data_server,
                                                exclude_anomalous=False,
                                                exclude_fmodel=False,
                                                keep_array_labels=True):
            if isinstance(label, list):
                if len(label)>3: continue
                all_labels.append( ",".join(ss.split(',')[0] for ss in label) )
            else:
                all_labels.append( label )


        if len(all_labels)<2:
            labin=all_labels[-1]
            self._logger( " ==> It seems that there is only one map available. Using %s" % labin )

        try:
            labin = labin.split(',')
            f = xray_data_server.get_xray_data(file_name = None, \
                                               labels = labin[0], \
                                               ignore_all_zeros = False, \
                                               parameter_scope="")

            phi = xray_data_server.get_xray_data(file_name = None, \
                                               labels = labin[1], \
                                               ignore_all_zeros = False, \
                                               parameter_scope="")

            if len(labin)>2:
                fom = xray_data_server.get_xray_data(file_name = None, \
                                               labels = labin[2], \
                                               ignore_all_zeros = False, \
                                               parameter_scope="")
            else:
                fom = None

            self.map_coeffs = iotbx.map_tools.combine_f_phi_and_fom(f=f, phi=phi, fom=fom)

        except:
            self._logger( " ==> Available MAP labels:" )
            for label in all_labels:
                if not isinstance(label, list): print( "         %s" % label )
            exit(0)

        self.map_coeffs.set_info(str(self.map_coeffs.info()).split('/')[-1])
        if verbose: self.map_coeffs.show_summary(prefix="     --> ")

        self.target_map = self.map_coeffs.fft_map(resolution_factor = 1.0/3.0)
        self.target_map.apply_sigma_scaling()


        self.target_miller_set = self.map_coeffs.miller_set(anomalous_flag = False,\
                                                            indices = self.map_coeffs.indices())

        self.map_data = self.target_map.real_map_unpadded()
        self.asu_mappings=self.symm.asu_mappings(buffer_thickness=0)

        assert reflection_file.file_content().n_symmetry_matrices() > 0

        s = maptbx.statistics(self.map_data)

        resol = self.target_miller_set.d_min()


        self.uc = self.symm.unit_cell()

        self._logger( " ==> Input map resolution is %.2f" % resol )

        self._logger( "     Space group: %s" % self.symm.space_group_info().symbol_and_number() )
        self._logger( "     Unit cell:   %s" % ",".join( [ "%.2f" % x for x in self.symm.unit_cell().parameters()] ) )
        self._logger( "     mean/sigma:  %.2f/%.2f" % (s.mean(),s.sigma()) )

        self._logger( "" )

        # the map reolution should be not higher than 2.0A
        if resol<2.0:
            self._logger( " ==> Filtering input map" )

            crystal_gridding              = maptbx.crystal_gridding(
                    unit_cell             = self.symm.unit_cell(),
                    space_group_info      = self.symm.space_group_info(),
                    symmetry_flags        = maptbx.use_space_group_symmetry,
                    pre_determined_n_real = self.target_map.n_real)

            complete_set = miller.build_set(
                            crystal_symmetry = self.symm,
                            anomalous_flag   = False,
                            d_min = 2.0)

            f_obs_cmpl = complete_set.structure_factors_from_map(
                                    map            = self.map_data,
                                    use_scale      = True,
                                    anomalous_flag = False,
                                    use_sg         = True)

            mtz_dataset = f_obs_cmpl.as_mtz_dataset(column_root_label="F")


            self.target_map = f_obs_cmpl.fft_map(resolution_factor = 1.0/3.0)
            self.target_map.apply_sigma_scaling()


            self.map_data = self.target_map.real_map_unpadded()



        if False:
            iotbx.ccp4_map.write_ccp4_map(
                file_name="random.map",
                unit_cell=self.symm.unit_cell(),
                space_group=self.symm.space_group_info().group(),
                gridding_first=self.map_data.accessor().origin(),
                gridding_last=self.map_data.accessor().last(),
                map_data=self.map_data,
                labels=flex.std_string(["iotbx.ccp4_map.tst"]))
            exit(1)


    # -------------------------------------------------------------------------


    def read_cloud(self, ifname=None):


        if ifname is None:
            ifname = self.template_fname

        with open(ifname, 'r') as ifile:
            pdb_string = ifile.read()

        ph = iotbx.pdb.input(source_info=None, lines=pdb_string).construct_hierarchy(sort_atoms=False)

        sel_cache = ph.atom_selection_cache()
        isel = sel_cache.iselection

        # get cloud
        sel_str = 'element H'
        ph_cloud = ph.select(isel(sel_str))


        # get backbone
        sel_str = "not element H"
        ph_bb = ph.select(isel(sel_str))

        return ph_bb, ph_cloud

    # -------------------------------------------------------------------------


    def local_region_density_correlation(self, large_unit_cell, \
                                        large_d_min, \
                                        large_density_map, \
                                        frag, \
                                        atom_radius=2.0):


        atoms = frag.atoms()
        atoms.set_b(new_b=flex.double(atoms.size(), 20.0))
        atoms.set_occ(new_occ=flex.double(atoms.size(), 1.0))

        sites_cart = frag.atoms().extract_xyz()

        sites_frac_large = large_unit_cell.fractionalize(sites_cart)

        large_n_real = large_density_map.focus()
        large_ucp = large_unit_cell.parameters()
        large_frac_min = sites_frac_large.min()
        large_frac_max = sites_frac_large.max()

        small_n_real = [0,0,0]
        small_origin_in_large_grid = [0,0,0]
        small_abc = [0,0,0]
        sites_frac_shift = [0,0,0]
        for i in range(3):
            grid_step = large_ucp[i] / large_n_real[i]
            buffer = 5.0
            grid_min = ifloor(large_frac_min[i] * large_n_real[i] - buffer)
            grid_max = iceil(large_frac_max[i] * large_n_real[i] + buffer)
            min_grid = grid_max - grid_min + 1
            small_n_real[i] = fftpack.adjust_gridding(min_grid=min_grid, max_prime=5)
            if (small_n_real[i] < large_n_real[i]):
                shift_min = (small_n_real[i] - min_grid) // 2
                small_origin_in_large_grid[i] = grid_min - shift_min
                small_abc[i] = small_n_real[i] * grid_step
                sites_frac_shift[i] = float(small_origin_in_large_grid[i]) / float(large_n_real[i])
            else:
                small_n_real[i] = large_n_real[i]
                small_origin_in_large_grid[i] = 0
                small_abc[i] = large_ucp[i]
                sites_frac_shift[i] = 0


        sites_cart_shift = large_unit_cell.orthogonalize(sites_frac_shift)
        sites_cart_small = sites_cart - sites_cart_shift


        small_xray_structure = xray.structure(
                                        crystal_symmetry=crystal.symmetry(
                                        unit_cell=tuple(small_abc)+large_ucp[3:], \
                                        space_group_symbol="P1"), \
                                        scatterers=frag.extract_xray_structure().scatterers())

        small_xray_structure.set_sites_cart(sites_cart=sites_cart_small)
        small_f_calc = small_xray_structure.structure_factors(
          d_min=large_d_min).f_calc()
        small_gridding = maptbx.crystal_gridding(unit_cell=small_f_calc.unit_cell(),
                                          space_group_info=small_f_calc.space_group_info(),
                                          pre_determined_n_real=small_n_real)

        small_fft_map = miller.fft_map(crystal_gridding=small_gridding,
                                        fourier_coefficients=small_f_calc)
        small_fft_map.apply_sigma_scaling()
        small_map = small_fft_map.real_map_unpadded()

        small_copy_from_large_map = maptbx.copy(
                        map_unit_cell=large_density_map,
                        first=small_origin_in_large_grid,
                        last=matrix.col(small_origin_in_large_grid) + matrix.col(small_n_real) - matrix.col((1,1,1)))

        assert small_copy_from_large_map.all() == small_map.all()




        result = []
        for rg in frag.residue_groups():

            grid_indices = maptbx.grid_indices_around_sites(
                        unit_cell=small_xray_structure.unit_cell(),
                        fft_m_real=small_n_real,
                        fft_n_real=small_n_real,
                        site_radii=flex.double(rg.atoms().size(), atom_radius),
                        sites_cart=flex.vec3_double(rg.atoms().extract_xyz()) - sites_cart_shift)

            corr = flex.linear_correlation( x=small_map.select(grid_indices), \
                                            y=small_copy_from_large_map.select(grid_indices))


            try:
                assert corr.is_well_defined()
                result.append(corr.coefficient())
            except:
                result.append()


        return result

    # -------------------------------------------------------------------------

    def density_score(self, sites_cart, frag=None):


        if frag is None:


            sites_frac = self.symm.unit_cell().fractionalize(sites_cart)

            score = []
            # do not refine backbone atoms!
            for site_frac in sites_frac[:-4]:
                score.append( self.map_data.eight_point_interpolation(site_frac) )

            return sum(score)

        else:

            rscc = self.local_region_density_correlation(large_unit_cell = self.symm.unit_cell(), \
                                                       large_density_map=self.map_data, \
                                                       large_d_min=self.dmin, \
                                                       frag = frag)[0]


            return rscc

    # -------------------------------------------------------------------------

    def local_fit_rotamer(self, residue, grm=None):


        def _refine(iresidue, imap_data, igrm, max_iterations=50):

            lbfgs_termination_params=scitbx.lbfgs.termination_parameters( max_iterations = max_iterations )

            minimized = real_space_refinement_simple.lbfgs(
                                density_map                 = imap_data,
                                sites_cart                  = iresidue.atoms().extract_xyz(),
                                gradients_method            = "fd",
                                real_space_target_weight    = (1.0),
                                lbfgs_termination_params    = lbfgs_termination_params,
                                geometry_restraints_manager = igrm,
                                real_space_gradients_delta  = 2.0/3.0)


            iresidue.atoms().set_xyz(new_xyz=minimized.sites_cart)


        # ---------------------------------------------------------------------


        tmp_frag = iotbx.pdb.hierarchy.root()
        tmp_frag.append_model(iotbx.pdb.hierarchy.model(id="0"))
        tmp_frag.models()[0].append_chain(iotbx.pdb.hierarchy.chain(id="0"))
        tmp_frag.only_chain().append_residue_group( residue.detached_copy() )

        sites_cart=residue.atoms().extract_xyz().deep_copy()
        resname = residue.only_atom_group().resname.strip().upper()
        comp_comp_id = self.mon_lib_srv.comp_comp_id_dict[resname.upper()]


        rotamer_iterator = comp_comp_id.rotamer_iterator(
                                    fine_sampling = True,
                                    atom_names    = residue.atoms().extract_name(),
                                    sites_cart    = sites_cart)


        if not rotamer_iterator.rotamer_info:
            _score = self.density_score(sites_cart, frag=tmp_frag)
            return _score


        hiscore = -666
        if rotamer_iterator.problem_message:
            print("WARNING: %s/%d: %s"%(residue.parent().id,residue.resseq_as_int(),rotamer_iterator.problem_message))
            return -3

        for r, rotamer_sites_cart in rotamer_iterator:
            tmp_frag.only_residue_group().atoms().set_xyz(new_xyz=rotamer_sites_cart)


            if grm:
                _refine(tmp_frag, self.map_data, grm, max_iterations=20)

            _score = self.density_score(rotamer_sites_cart, frag=tmp_frag)

            if _score>hiscore:
                hiscore = _score
                best_rotamer_sites_cart = rotamer_sites_cart

        residue.atoms().set_xyz(new_xyz=best_rotamer_sites_cart)


        # final, a bit longer refinement (better leave it to a more sphisticated methods)
        #_refine(residue, self.map_data, grm, max_iterations=50)


        # calc RSCC
        return self.density_score(residue.atoms().extract_xyz(), frag=tmp_frag)

    # -------------------------------------------------------------------------


    def generate_restraints_manager_for_ag(self, ag, verbose=False):

        """
            a simple(r) replacement to cctbx restraints manager based on geostd library
        """


        def i_atom_by_id(name):
            return atom_dict[name]



        atom_dict = dict( [ (atom.name.strip().upper(), idx) for idx,atom in enumerate(ag.atoms()) ] )



        # restraints include hydrogens
        comp = self.mon_lib_srv.comp_comp_id_dict[ag.resname.upper()]
        motif = comp.as_geometry_restraints_motif()


        bond_proxies = geometry_restraints.bond_sorted_asu_proxies(asu_mappings=None)
        angle_proxies = geometry_restraints.shared_angle_proxy()


        dihedral_proxies = geometry_restraints.shared_dihedral_proxy()
        chirality_proxies = geometry_restraints.shared_chirality_proxy()
        planarity_proxies = geometry_restraints.shared_planarity_proxy()



        # add covalent bonds
        for bond in motif.bonds_as_list():

            # again, rough solution for the missing H problem
            try:
                i_seqs = tuple(map(i_atom_by_id, bond.atom_names))
            except:
                continue

            bond_proxies.process(geometry_restraints.bond_simple_proxy(
                                        i_seqs         = i_seqs, \
                                        distance_ideal = bond.distance_ideal, \
                                        weight         = bond.weight))

        # add angles
        for angle in motif.angles_as_list():

            # rough solution for the missing H problem
            try:
                i_seqs = tuple(map(i_atom_by_id, angle.atom_names))
            except:
                continue

            angle_proxies.append(geometry_restraints.angle_proxy(
                                        i_seqs      = i_seqs, \
                                        angle_ideal = angle.angle_ideal, \
                                        weight      = angle.weight))

        # add planes (if any)
        # this one is special, since plane defs may contain lots of hydrogens, and still be valid
        for plane in motif.planarities_as_list():

            weights = []
            i_seqs = []
            for _w, _n in zip( plane.weights, plane.atom_names ):

                try:
                    i_seqs.append( i_atom_by_id(_n) )
                except:
                    continue

                weights.append(_w)

            planarity_proxies.append(geometry_restraints.planarity_proxy(
                                        i_seqs         = i_seqs, \
                                        weights        = weights))




        # add chiralities
        for chir in motif.chiralities_as_list():
            chirality_proxies.append(geometry_restraints.chirality_proxy(
                                        i_seqs         = tuple(map(i_atom_by_id, chir.atom_names)), \
                                        both_signs     = chir.both_signs,   \
                                        volume_ideal   = chir.volume_ideal, \
                                        weight         = chir.weight))


        # add dihedrals
        for dih in motif.dihedrals_as_list():
            try:
                i_seqs = tuple(map(i_atom_by_id, dih.atom_names))
            except:
                continue


            dihedral_proxies.append(geometry_restraints.dihedral_proxy(
                                        i_seqs         = i_seqs,          \
                                        angle_ideal    = dih.angle_ideal, \
                                        periodicity    = dih.periodicity, \
                                        weight         = dih.weight))



        bond_params_table = geometry_restraints.extract_bond_params(
                                n_seq=ag.atoms().size(),
                                bond_simple_proxies=bond_proxies.simple)


        grmanager = geometry_restraints.manager.manager(
                                        crystal_symmetry = self.symm,           \
                                        bond_params_table = bond_params_table,  \
                                        angle_proxies     = angle_proxies,      \
                                        planarity_proxies = planarity_proxies,  \
                                        chirality_proxies = chirality_proxies,  \
                                        dihedral_proxies  = dihedral_proxies)



        if verbose:
            msg = " ==> [%s] : %3i angle, %3i bond, %3i dihedral, and %3i planarity proxies" % \
                                                                (ag.resname,\
                                                                 grmanager.get_n_angle_proxies(),\
                                                                 grmanager.get_n_bond_proxies(),\
                                                                 grmanager.get_dihedral_proxies().size(),\
                                                                 grmanager.get_n_planarity_proxies())



        return grmanager


    # -------------------------------------------------------------------------

    def mutate_aa(self, res, target_aa, target_resseq=None, refine=True, verbose=False):


        old_res = res

        # build ALA instead of UNKS
        if target_aa in ['UNK']:
            ideal_ph = self.ideal_dict['ala'].deep_copy()
        else:
            ideal_ph = self.ideal_dict[target_aa.lower()].deep_copy()

        ideal_ag = ideal_ph.only_model().only_chain().only_residue_group().only_atom_group()

        tmp_residue = old_res.only_atom_group().detached_copy()

        ideal_xyz = flex.vec3_double(self.res_bb(ideal_ag))
        try:
            old_xyz = flex.vec3_double(self.res_bb(old_res))
        except:
            return None

        superposition = superpose.least_squares_fit(old_xyz, ideal_xyz, method=["kearsley", "kabsch"][0])

        rtmx = matrix.rt((superposition.r, superposition.t))
        ideal_ag.atoms().set_xyz( rtmx * ideal_ag.atoms().extract_xyz() )


        old_res_ag_copy = old_res.only_atom_group().detached_copy()
        old_res.remove_atom_group(old_res.only_atom_group())
        old_res.append_atom_group(ideal_ag.detached_copy())


        # preserve backbone atom coordinates
        for atom in old_res.atoms():
            atom.set_b(30)
            if(atom.name.strip().upper() in ["N", "C", "CA", "O"]):
                old_res.only_atom_group().remove_atom(atom)

        for atom in old_res_ag_copy.atoms():
            if(atom.name.strip().upper() in ["N", "C", "CA", "O"]):
               old_res.only_atom_group().append_atom(atom.detached_copy())

        # do not refine UNKs
        if target_aa.upper() in ['UNK']:
            return 0

        # refine side-chain atoms
        if refine:
            ideal_ph.only_residue_group().remove_atom_group( ideal_ph.only_atom_group() )
            ideal_ph.only_residue_group().append_atom_group( old_res.only_atom_group().detached_copy() )


            _restraints_manager = self.generate_restraints_manager_for_ag( old_res.only_atom_group(), verbose=False )

            reference_sites = flex.vec3_double()
            reference_selection = flex.size_t()
            cntr = 0
            for atom in old_res.atoms():
                if atom.name.strip().upper() in ["N", "C", "CA", "O"]:
                    reference_sites.append(atom.xyz)
                    reference_selection.append(cntr)
                cntr += 1



            _restraints_manager.adopt_reference_coordinate_restraints_in_place( reference.add_coordinate_restraints( sites_cart = reference_sites, \
                                                                                                                      selection = reference_selection, \
                                                                                                                          sigma = 0.01) )

        else:
            _restraints_manager=None

        _r = self.local_fit_rotamer(old_res, _restraints_manager)

        if target_resseq: old_res.resseq = target_resseq%10000

        return max(0,_r)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
