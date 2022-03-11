#!/usr/bin/env python

'''
Secondary structure
Repository
https://github.com/cmbi/dssp

Install dependencies
(Ubuntu) sudo apt-get install libboost-all-dev autoconf
(MacOS) Install with Homebrew: autoconf, automake, boost

Compile:
./autogen.sh
./configure
make mkdssp
./mkdssp (to test)
'''

import math
import numpy as np
from Bio.PDB import PDBList, DSSP, HSExposureCB, PPBuilder
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
from matplotlib import patches
from Bio.SeqUtils import seq1


rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                  (-180, 50, 80, 130, 'E', 'blue'),
                  (-100, -180, 100, 60, 'P', 'green'),
                  (-100, 50, 100, 130, 'P', 'green'),
                  (-180, -120, 180, 170, 'H', 'red'),
                  (0, -180, 180, 360, 'L', 'yellow')]


# Ramachandran regions
regions_matrix = []
with open("ramachandran.dat") as f:
    for line in f:
        if line:
            regions_matrix.append([int(ele) for ele in line.strip().split()])


if __name__ == "__main__":

    pdb_ids = [('2kkw', 'A')]
    for pdb_id, chain_id in pdb_ids:

        # Fetch a PDB file to the current dir
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(pdb_id, pdir='data/', file_format='pdb')  # Will save to pdbXXXX.ent

        # Load the structure
        structure = PDBParser(QUIET=True).get_structure(pdb_id, "data/pdb{}.ent".format(pdb_id))

        data = []
        for model in structure:
            print(model[chain_id])

            model_data = []

            residues = [residue for residue in model[chain_id]]

            dssp = DSSP(model, "data/pdb{}.ent".format(pdb_id), dssp="binx/dssp/mkdssp")  # WARNING Check the path of mkdssp
            dssp = dict(dssp)
            hse = HSExposureCB(model[chain_id])

            ppb = PPBuilder()  # PolyPeptideBuilder
            rama_dict = {}
            for pp in ppb.build_peptides(model[chain_id]):
                phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
                for i, residue in enumerate(pp):
                    phi, psi = phi_psi[i]
                    ss_class = None
                    if phi is not None and psi is not None:
                        for x, y, width, height, ss_c, color in rama_ss_ranges:
                            if x <= phi < x + width and y <= psi < y + height:
                                ss_class = ss_c
                                break
                    rama_dict[(chain_id, residue.id)] = [phi, psi, ss_class]


            for residue in residues:
                model_data.append((residue, dssp.get((chain_id, residue.id))[2],
                                   dssp.get((chain_id, residue.id))[3],
                                    *hse[(chain_id, residue.id)],
                                   *rama_dict.get((chain_id, residue.id))))

            # Transpose elements
            data.append(list(map(list, zip(*model_data))))

        # Transpose elements
        data = list(map(list, zip(*data)))

        fig, axes = plt.subplots(5, 1, figsize=(16, 24))

        c = 0
        for i, (feature, name) in enumerate(zip(data, ['residue', 'ss', 'asa', 'hse1', 'hse2', 'hse3', 'phi', 'psi', 'ss_three'])):
            if i in [2, 6, 7]:

                feature = np.array(feature, dtype=np.float)
                mean = np.nanmean(feature, axis=0)
                std = np.nanstd(feature, axis=0)

                # print(mean)
                # print(std)

                axes[c].errorbar(np.arange(len(mean)), mean, yerr=std, ecolor='r', capthick=0.5)

                axes[c].set_title(name)
                c += 1

        plt.tight_layout()  # Remove figure padding
        plt.savefig('data/{}_features.png'.format(pdb_id), bbox_inches='tight')