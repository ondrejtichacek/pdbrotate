#!/bin/python

import os
import sys
import getopt

import numpy as np
from scipy.linalg import expm3, norm

import mdtraj as md

HELP = """
DESCRIPTION

Rotates molecule around axis.

OPTIONS

        -f [.pdb]  input file
        -o [.pdb]  output file

   --angle []      angle of rotation in degrees (float number)
    --axis []      axis of rotation (comma separated list of numbers)

        -v         verbose flag (--verbose)
        -h         print this help (--help)

EXAMPLE USAGE

Rotate the molecule 90 degrees around the z axis.

python3 rotate.py -f peptide.pdb -o out.pdb --angle=90 --axis=0,0,1
"""

def usage():
    print(HELP)

def parse_cmd_options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:o:v:", ["help", "verbose", "angle=", "axis="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    foutput = None
    finput = None
    verbose = False
    angle = None
    axis = None

    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-f"):
            finput = a
        elif o in ("-o"):
            foutput = a
        elif o in ("--angle"):
            angle = [float(a),]
        elif o in ("--axis"):
            axis = [np.asarray([float(x) for x in a.split(",")]),]
        else:
            assert False, "unhandled option"

    theta = [ang * np.pi / 180 for ang in angle]

    return {
        'theta': theta,
        'axis': axis,
        'verbose': verbose,
        'finput': finput,
        'foutput': foutput,
    }

def rotation_matrix(axis, theta):
    """http://stackoverflow.com/questions/6802577/"""
    return expm3(np.cross(np.eye(3), axis/norm(axis)*theta))

def move(traj, fun):
    newcoordinates = []
    for v in traj.xyz[0]:
        newcoordinates.append(fun(v))
    traj.xyz = [newcoordinates]

def rotate(finput, foutput, theta, axis, verbose):

    top = md.load_topology(finput)
    atoms_to_load = top.select("all")
    atoms = md.load(
        finput,
        atom_indices=atoms_to_load)

    coordinates = atoms.xyz[0]

    center_of_mass = md.compute_center_of_mass(atoms)[0]

    for th, ax in zip(theta, axis):

        M = rotation_matrix(ax, th)
        move(atoms, lambda v: np.dot(M, v))

    new_center_of_mass = md.compute_center_of_mass(atoms)[0]

    move(atoms, lambda v: v - new_center_of_mass + center_of_mass)

    atoms.save(foutput)


if __name__ == "__main__":

    rotate(**parse_cmd_options())
