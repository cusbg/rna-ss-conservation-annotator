# rna-ss-conservation-annotator
## A tool for secondary structure-based annotation of a set of RNA molecules


Given multiple sequences with its secondary structures,
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

Run:
  >>> compare.py --help
  
to get a description of all command-line parameters.

