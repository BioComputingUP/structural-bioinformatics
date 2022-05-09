#!/usr/bin/env python

# **** WARNING ****
# For some Ubuntu versions is not possible to install Pymol with Conda

# To run this script with the Pymol version installed on your system you have to exit the virtual environment
# for Conda just type "conda deactivate"

# **** WARNING (MacOS only) ****
# PyMOL and Python scripts only works if installed with the Homebrew package manager

# To execute a Python script (like this one):
# python pymol_lines.py


import pymol
from pymol import cgo, cmd, util

pymol.finish_launching()  # Open Pymol

# Input 1jsu (3 chains, non-standard amino acids), 1az5 (disordered loops, chain breaks)
pdb_id = '1jsu'

cmd.fetch(pdb_id, pdb_id, path="data/")  # Download the PDB
# cmd.load("data/pdb{}.ent".format(pdb_id), pdb_id)  # Load from file

cmd.remove("resn hoh")  # Remove water molecules
cmd.hide("lines", "all")  # Hide lines
cmd.show("cartoon", pdb_id)  # Show cartoon
cmd.show("sticks", "hetatm")  # Show hetero atoms as sticks
# # cmd.spectrum(selection="all")  # Rainbow color
util.cbc(selection="all")  # Color by chain
util.cnc(selection="all")  # Color by atom type, but not the C atoms

# Select and color two residues
sele_name = "nodes"
res1 = 'B/200/'
res2 = 'C/52/'
cmd.select(sele_name, '{} or {}'.format(res1, res2))
cmd.show("spheres", sele_name)
cmd.set('sphere_transparency', 0.5, sele_name)  # Set transparency

# Get coordinates of two atoms
atom1 = 'B/200/SD'
atom2 = 'C/52/CE'
coord1 = cmd.get_coords(atom1, 1)  # shape N * 3, where N is the number of atoms in the selection
coord2 = cmd.get_coords(atom2, 1)
# model = cmd.get_model(pdb_id, 1)  # Fastest way to get all atom coordinates

# Calculate center of mass between two residues and create a new "pseudo" atom
center_of_mass = (coord1[0] + coord2[0]) / 2
print(coord1, coord2, center_of_mass)

obj_name = "ps_atom"
cmd.pseudoatom(obj_name, pos=list(center_of_mass))
cmd.show("spheres", obj_name)
# cmd.extract(...  # Move selected atoms to a new object
# cmd.create(...  # Create a new object from selection

cr, cg, cb = (1.0, 0.0, 0.0)  # RGB red

# Create lines object
obj = [cgo.BEGIN, cgo.LINES, cgo.COLOR, cr, cg, cb]
obj.append(cgo.VERTEX)
obj.extend(list(coord1[0]))
obj.append(cgo.VERTEX)
obj.extend(list(coord2[0]))
obj.append(cgo.END)

# Set the object
obj_name = 'edges'
cmd.load_cgo(obj, obj_name)
cmd.set("cgo_line_width", float(3), obj_name)

cmd.orient(pdb_id)  # Set the origin to full protein


# *** Exercises ****

# *** Identify chain breaks in 1az5
# *** Draw a line connecting broken chain fragments
# *** Draw a line connecting the center of mass of each residue (assume atom mass equal to 1)
