#! /usr/bin/env python
"""
condlib.py
Functions for generating conductance matrices for crossbar memory arrays

Created by Jeremy Smith on 2017-07-06
University of California, Berkeley
"""

import numpy as np

__author__ = "Jeremy Smith"
__version__ = "1.0"


def conductance_matrix_WRITE(n, rline, rcell, vWLsel, vWLnsel, vBLsel, vBLnsel,
                             isel=0, jsel=0, verbose=True):
    """
    Calculates the conductance matrix and the vector of supply currents for a
    WRITE event i.e. all word and bit lines have an applied voltage
        n x n matrix
        power supplied on 0-index side of array (single side supply)
        word line node listed followed by bit line nodes
        no pull up resistor (equal to line resistance on power line side)
    """

    # Conductance coefficients
    a2 = 1.0 / rcell + 2.0 / rline
    a1 = 1.0 / rcell + 1.0 / rline
    b = -1.0 / rline
    c = -1.0 / rcell
    ds = vWLsel / rline
    dn = vWLnsel / rline

    # Setup empty conductance matrix and current vector
    nsq = n * n
    mat = np.zeros((2 * nsq, 2 * nsq))
    iin = np.zeros(2 * nsq)

    # Loop for word line nodes (W) and bit line nodes (B)
    for l in ['W', 'B']:
        # Loop for all word lines (i)
        for i in xrange(n):
            # Loop for all bit lines (j)
            for j in xrange(n):
                if l is 'W':
                    p = n*i + j
                    if verbose:
                        print l, i, j, p
                    if j == 0:           # i.e. power line side
                        mat[p][p] = a2
                        mat[p][p + 1] = b
                        mat[p][p + nsq] = c
                        iin[p] = ds if i is isel else dn

                    elif j == (n - 1):   # i.e. floating side
                        mat[p][p] = a1
                        mat[p][p - 1] = b
                        mat[p][p + nsq] = c

                    else:                # i.e. middle of array
                        mat[p][p] = a2
                        mat[p][p - 1] = b
                        mat[p][p + 1] = b
                        mat[p][p + nsq] = c

                else:
                    p = n*i + j + nsq
                    if verbose:
                        print l, i, j, p
                    if i == 0:           # i.e. power line side
                        mat[p][p] = a2
                        mat[p][p + n] = b
                        mat[p][p - nsq] = c
                        iin[p] = ds if j is jsel else dn

                    elif i == (n - 1):   # i.e. floating side
                        mat[p][p] = a1
                        mat[p][p - n] = b
                        mat[p][p - nsq] = c

                    else:                # i.e. middle of array
                        mat[p][p] = a2
                        mat[p][p - n] = b
                        mat[p][p + n] = b
                        mat[p][p - nsq] = c

    return mat, iin


def conductance_matrix_READ(n, rline, rcell, rpu, vWLsel, vBLsel,
                            isel=0, jsel=0, verbose=True):
    """
    Calculates the conductance matrix and the vector of supply currents for a
    standard READ event i.e. selected word and bit lines have an applied
    voltage and others are floating
        n x n matrix
        power supplied on 0-index side of array (single side supply)
        word line node listed followed by bit line nodes
        pull up resistor (rpu) can be set for selected line
    """

    # Conductance coefficients
    a2 = 1.0 / rcell + 2.0 / rline
    a1 = 1.0 / rcell + 1.0 / rline
    a3 = 1.0 / rcell + 1.0 / rline + 1.0 / rpu
    b = -1.0 / rline
    c = -1.0 / rcell
    ds = vWLsel / rpu

    # Setup empty conductance matrix and current vector
    nsq = n * n
    mat = np.zeros((2 * nsq, 2 * nsq))
    iin = np.zeros(2 * nsq)

    # Loop for word line nodes (W) and bit line nodes (B)
    for l in ['W', 'B']:
        # Loop for all word lines (i)
        for i in xrange(n):
            # Loop for all bit lines (j)
            for j in xrange(n):
                if l is 'W':
                    p = n*i + j
                    if verbose:
                        print l, i, j, p
                    if j == 0:           # i.e. power line side
                        if i is isel:
                            mat[p][p] = a3
                            mat[p][p + 1] = b
                            mat[p][p + nsq] = c
                            iin[p] = ds
                        else:
                            mat[p][p] = a1
                            mat[p][p + 1] = b
                            mat[p][p + nsq] = c

                    elif j == (n - 1):   # i.e. floating side
                        mat[p][p] = a1
                        mat[p][p - 1] = b
                        mat[p][p + nsq] = c

                    else:                # i.e. middle of array
                        mat[p][p] = a2
                        mat[p][p - 1] = b
                        mat[p][p + 1] = b
                        mat[p][p + nsq] = c

                else:
                    p = n*i + j + nsq
                    if verbose:
                        print l, i, j, p
                    if i == 0:           # i.e. power line side
                        if j is jsel:
                            mat[p][p] = a3
                            mat[p][p + n] = b
                            mat[p][p - nsq] = c
                            iin[p] = ds
                        else:
                            mat[p][p] = a1
                            mat[p][p + n] = b
                            mat[p][p - nsq] = c

                    elif i == (n - 1):   # i.e. floating side
                        mat[p][p] = a1
                        mat[p][p - n] = b
                        mat[p][p - nsq] = c

                    else:                # i.e. middle of array
                        mat[p][p] = a2
                        mat[p][p - n] = b
                        mat[p][p + n] = b
                        mat[p][p - nsq] = c

    return mat, iin
