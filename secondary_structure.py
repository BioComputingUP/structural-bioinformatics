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


pdb_id = '1az5'

# Fetch a PDB file to the current dir
pdbl = PDBList()
pdbl.retrieve_pdb_file(pdb_id, pdir='data/', file_format='pdb') # Will save to pdbXXXX.ent

# Load the structure
structure = PDBParser(QUIET=True).get_structure(pdb_id, "data/pdb{}.ent".format(pdb_id))

# Secondary structure, solvent accesibility (with DSSP)
dssp = DSSP(structure[0], "data/pdb{}.ent".format(pdb_id), dssp="binx/dssp/mkdssp")  # WARNING Check the path of mkdssp
# for ss in dssp:
#     # dssp index, amino acid, secondary structure, relative ASA, phi, psi,
#     # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
#     # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy
#     print(ss)
#
dssp_dict = dict(dssp)
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             print(chain.id, residue.id, dssp_dict.get((chain.id, residue.id)))
#
# # HSE exposure
# # https://en.wikipedia.org/wiki/Half_sphere_exposure
# # Low HSE-up means high ASA (absolute solvent accessibility)
# hse = HSExposureCB(structure[0])
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             if residue.id[0] == " ":
#                 print(chain.id, residue.id, hse[(chain.id, residue.id)], dssp_dict.get((chain.id, residue.id)))  # HSE beta up (EXP_HSE_B_U), HSE beta down (EXP_HSE_B_D)


# Secondary structure from phi & psi
# https://pubs.acs.org/doi/10.1021/ja306905s#

# List of Ramachandran areas corresponding to different secondary structure classes
# E = beta sheet, P = polyproline I && II,
# H = alpha-helix, R = left-handed helix
# (lower-left phi and psi, width, height, class, color)
rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                  (-180, 50, 80, 130, 'E', 'blue'),
                  (-100, -180, 100, 60, 'P', 'green'),
                  (-100, 50, 100, 130, 'P', 'green'),
                  (-180, -120, 180, 170, 'H', 'red'),
                  (0, -180, 180, 360, 'L', 'yellow')]

# Calculate PSI and PHI
ppb = PPBuilder()  # PolyPeptideBuilder
rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }
for model in structure:
    for chain in model:
        for pp in ppb.build_peptides(chain):

            phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
            for i, residue in enumerate(pp):
                # print(model, chain, i, residue, phi_psi[i])

                # Convert radians to degrees and remove first and last value that are None
                if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                    rama.setdefault(chain.id, [[], [], []])
                    rama[chain.id][0].append(residue)
                    rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
                    rama[chain.id][2].append(math.degrees(phi_psi[i][1]))

# Get SS from phi/psi and compare with DSSP
for chain_id in rama:
    for residue, phi, psi in zip(*rama[chain_id]):
        ss_class = None
        for x, y, width, height, ss_c, color in rama_ss_ranges:
            if x <= phi < x + width and y <= psi < y + height:
                ss_class = ss_c
                break
        print(residue, ss_class, dssp_dict.get((chain_id, residue.id))[2], phi, psi)

# Plot Ramachandran SS regions
f, axes = plt.subplots(1, len(rama), figsize=(12, 12))
axes = np.array(axes).reshape(-1)  # Hack to give consistency for single/multiple suplots (-1 force reshape to infer dimensions)
for ax, chain_id in zip(axes, rama):

    # Plot SS regions
    for x, y, width, height, _, color in rama_ss_ranges:
        ax.add_patch(patches.Rectangle(
            (x, y),  # (x,y)
            width,  # width
            height,  # height
            color=color, zorder=0))  # transparency

    # Plot points
    ax.scatter(rama[chain_id][1], rama[chain_id][2], s=6, color='black', zorder=1)

    ax.set_xlabel('phi')
    ax.set_ylabel('psi')

plt.tight_layout()  # Remove figure padding
plt.savefig('data/{}_ramachandran_ss.png'.format(pdb_id), bbox_inches='tight')
