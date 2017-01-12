#!/bin/python

import os
import sys
import getopt

import numpy as np
from scipy.linalg import expm3, norm

import mdtraj as md

from rotate import rotate

HELP = """
DESCRIPTION

Aligns the molecule according to its principial axis.

OPTIONS

        -f [.pdb]  input file
        -o [.pdb]  output file

        -v         verbose flag (--verbose)
        -h         print this help (--help)

EXAMPLE USAGE

python3 align.py -f peptide.pdb -o out.pdb
"""

def usage():
    print(HELP)

def parse_cmd_options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:o:v:", ["help", "verbose"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    foutput = None
    finput = None
    verbose = False

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
        else:
            assert False, "unhandled option"

    return {
        'verbose': verbose,
        'finput': finput,
        'foutput': foutput,
    }

from math import acos, atan2, cos, sin
from numpy import array, float64, zeros
from numpy.linalg import norm


def cartesian_to_spherical(vector):
    """Convert the Cartesian vector [x, y, z] to spherical coordinates [r, theta, phi].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.


    @param vector:  The Cartesian vector [x, y, z].
    @type vector:   numpy rank-1, 3D array
    @return:        The spherical coordinate vector [r, theta, phi].
    @rtype:         numpy rank-1, 3D array
    """

    # The radial distance.
    r = norm(vector)

    # Unit vector.
    unit = vector / r

    # The polar angle.
    theta = acos(unit[2])

    # The azimuth.
    phi = atan2(unit[1], unit[0])

    # Return the spherical coordinate vector.
    return array([r, theta, phi], float64)


def spherical_to_cartesian(spherical_vect, cart_vect):
    """Convert the spherical coordinate vector [r, theta, phi] to the Cartesian vector [x, y, z].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.


    @param spherical_vect:  The spherical coordinate vector [r, theta, phi].
    @type spherical_vect:   3D array or list
    @param cart_vect:       The Cartesian vector [x, y, z].
    @type cart_vect:        3D array or list
    """

    # Trig alias.
    sin_theta = sin(spherical_vect[1])

    # The vector.
    cart_vect[0] = spherical_vect[0] * cos(spherical_vect[2]) * sin_theta
    cart_vect[1] = spherical_vect[0] * sin(spherical_vect[2]) * sin_theta
    cart_vect[2] = spherical_vect[0] * cos(spherical_vect[1])

def rotation_matrix_2d(theta):

    return [[cos(theta), -sin(theta)],
            [sin(theta),  cos(theta)]]

    # return [[cos(theta), -sin(theta), 0],
    #         [sin(theta),  cos(theta), 0],
    #         [0,           0,          1]]

def align(finput, foutput, verbose):

    top = md.load_topology(finput)
    atoms_to_load = top.select("name CA")
    atoms = md.load(
        finput,
        atom_indices=atoms_to_load)

    # vec = atoms.xyz[0][-1] - atoms.xyz[0][0]

    data = atoms.xyz[0]

    # Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean = data.mean(axis=0)

    # Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean)

    vec = vv[0]

    print("PRINCIPIAL AXIS: {}".format(vec))

    spherical = cartesian_to_spherical(vec)
    theta = spherical[1]
    phi = spherical[2]

    print("THETA: {}".format(theta*180/np.pi))
    print("PHI: {}".format(phi*180/np.pi))

    rotate(
        axis=[[0,0,1], [0,1,0]],
        theta=[-phi, np.pi - theta],
        finput=finput,
        foutput=foutput,
        verbose=verbose)

if __name__ == "__main__":

    align(**parse_cmd_options())
