#!/usr/bin/env python

"""
Comparing secondary structure conservancy
======================

Given:

* multiple sequences with its secondary structures,

produces a graphical representation of secondary structures conservation
comparation with respect to one of the sequences.

Expects directory containing fasta files (.fst) which includes (predicted)
secondary structure and a list file, that contains names of the sequences
that should be compared. First sequence in the list is always the primary,
meaning all other sequences are compared to it.


Requirements
------------

Requires the following to be installed:

* The ``rPredictor`` package libraries (which should have been installed
  using the rPredictor package ``setup.py``)

* ``RNAPlot``, which can be downloaded as a part of the `Vienna RNA package <http://www.tbi.univie.ac.at/RNA/>`_

* ``clustalw2``, which can be found somewhat unorthodoxly at the `Help page of its website <http://www.ebi.ac.uk/Tools/msa/clustalw2/help/>`_


Usage
-----

The ``compare.py`` script can be used as a standalone command-line tool.

Command-line usage
^^^^^^^^^^^^^^^^^^

Run::

  >>> compare.py --help

to get a description of all command-line parameters.

-------------------
"""

import logging
import argparse
import StructureComparator
import os


def _build_argument_parser():
    """Creates the argument parser for running the :func:`main` function.

    :rytpe: argparse.ArgumentParser
    :returns: The argument parser to which to pass the command-line options for
        ``compare.py``.
    """

    parser = argparse.ArgumentParser(description='Comparation of conservancy le' +
                                     'vels among multiple sequences.',
                                     add_help=True)

    parser.add_argument("sequences_dir", metavar='SEQ_DIR',
                        help='Path to directory containing fasta files with seq' +
                        'uences to be compared. The directory should also conta' +
                        'in list file specifying the sequence files.')
    parser.add_argument("-l", "--list_file", metavar='LIST_FILE', default='list.txt',
                        help='File containing names of sequences that should be' +
                        'compared. All listed sequences must be placed in a fil' +
                        'e with the same name in a same dir.')
    parser.add_argument("-w", "--winsizes", metavar='WINSIZES_LIST',
                        default='20,40,60',
                        help='List of comma-separated numbers giving which wind' +
                        'ow sizes should be used for comparison.')
    parser.add_argument("-e", "--segments_file", action='store', default=None,
                        help='Name of the file containing expansion segments of' +
                        ' some sequences from the set. File must contain valid ' +
                        'Python dictionary indexed by names of the sequences (a' +
                        's in list file).')
    parser.add_argument("-s", "--stat_file", action='store', default=None,
                        help='Name of the csv file for storing comparation stat' +
                        'istics. If not set, statistics are not stored.')

    return parser

def main(args):
    # if statistisc file name is given, prepare the output csv file
    if args.stat_file != None:
        statFileName = os.path.join(args.sequences_dir, args.stat_file)
        statFile = open(statFileName, 'w')
        statFile.write("Name,Win. size,Total mismatch,ES mismatch,Total vs ES ratio,Avg.mismatch,Avg. mis. in ES,Avg. mis. outside ES\n")
        statFile.close()
    else:
        statFileName = None
    # loop over different windows and generate comparations for each window size
    winsize_list = [int(x) for x in args.winsizes.split(',') if x.strip().isdigit()]
    for winSize in winsize_list:
        # prepare comparator
        comparator = StructureComparator.ConservancyComparator(
          args.sequences_dir, args.list_file, winSize, 
          statFileName, args.segments_file)
        # run comparaison
        comparator.compare()
        print "Finished comparing for window size: " + str(winSize)

if __name__ == '__main__':
    parser = _build_argument_parser()
    args = parser.parse_args()
    main(args)
