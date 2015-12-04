#
# StructureComparator.py
#
# Module for comparing multiple sequence structures rating their conservancy level against each other
#
"""This module provides ConservancyComparator class which compares
secondary structures conservancy to each other.

The sequences are loaded from a list file - that is a file containing names
of all sequences to compare, one name per line. The names in this list and
names of the fasta files containing the sequence must be the same. They also
need to be in the same directory (path to the directory is a parameter).

------------------

"""
__author__ = "Jan Pesek"
__contributors__ = ""
__credits__ = ["Jan Hajic", "David Hoksza", "Josef Panek"]
__maintainer__ = "Jan Pesek"
__email__ = "jpesek89@gmail.com"
__status__ = "Prototype"

from Bio.Align.Applications import ClustalwCommandline
from AlignmentTools import extract_seq_string_from_aln
from AlignmentTools import read_alignment
from ComparationTagger import ComparationTagger
import os
import re
import subprocess
import shutil
import shlex
import sys

class ConservancyComparator(object):
    """A ConservancyComparator is a class that implements the comparing
     and raning algorhitm for set of sequences with specific window size.

    Usage:

    >>> comparator = ConservancyComparator('sequences_dir', 'list.txt', 40, 0.5)
    >>> comparator.compare()

    etc.
    """

    def __init__(self, dir_path, list_name, window_size, statfile_name = None, segments_file = None):
        """Initializes the ConservancyComparator.

        :type dir_path: string
        :param dir_path: The directory which contains fasta files with sequences
            that should be compared (and also list file and possibly segments file).
            This directory will also contain class output (mainly svg files).

        :type list_name: string
        :param list_name: Name of the file that contains names of the sequences
            to be compared. Each line should be one filename, first record
            would be considered a primary sequence (the one to which the others
            are being compared to).

        :type window_size: int
        :param window_size: The number of consecutive positions
            that will be considered for calculating mismatch percentage.

        :type statfile_name: string
        :param statfile_name: Name of the output statistics file, if desired.
            The file will be placed in the dir_parh directory. This is especially
            useful in combination with segments_file parameter. Default is None.

        :type segments_file: string
        :param segments_file: Name of the file containing so called expansion segments.
            These are segments that should be outputted differently in the
            visualization. Compare method also calculates special statistics for
            these segments. The file should contain expansion segments as a list
            of Python Tuples in a Python Dictionary indexed by names
            of the sequences (as listed in list_name file). Default is None.
        """
        self.dir_path = dir_path
        self.window_size = window_size
        self.list_name = list_name
        self.mismatch_percentage = 0.5 # splits structure to regions by this value - can be used by other rPredictor components
        self.segments_file = segments_file
        self.statfile_name = statfile_name
        self.list_filename = os.path.join(self.dir_path, self.list_name)
        alignment_filename_helper = os.path.splitext(self.list_filename)[0]
        self.all_sequences_filename = alignment_filename_helper + '_all_sequences.txt'
        self.alignment_filename = alignment_filename_helper + '.aln'
        self.gap_open = 20
        self.gap_extension = 0.1

    def is_sequence(self, seq_string):
        """Checks, whether string is a valid sequence (NOTE: all 5 letters
        A, C, G, T, U are considered as a valid sequence).

        :type seq_string: string
        :param seq_string: String to be tested

        :rtype: bool
        :returns: True, if provided string is a sequence.
        """
        return bool(re.match("^[ACGTU]*$", seq_string))

    def _extract_converted_ss(self, input_handle):
        """Extracts secondary structure in a dot-paren notation from provided
            fasta file and converts it to a ClustalW-valid format.

        :type input_handle: FileObject
        :param input_handle: Fasta file with secondary structure

        :rtype: string
        :returns: Converted secondary structure (dot-paren to ClustalW ACG format)
        """
        header = input_handle.next()
        ln = input_handle.next()
        while(self.is_sequence(ln)):
            ln = input_handle.next()
        dp = ln
        try:
            while not self.is_sequence(ln):
                ln = input_handle.next()
                dp = dp + ln
        except StopIteration:
            pass
        return header + dp.replace('(', 'A').replace(')', 'C').replace('.', 'G') + '\n'

    def _create_matrix_file(self, filepath):
        """Creates simple file containing matrix to be used by ClustalW for aligning.

        :type filepath: string
        :param filepath: name (path) of the file
        """
        f = open(filepath, 'w')
        f.write("    A   G   C   U   *\n")
        f.write("A   1   0   0   0   0\n")
        f.write("G   0   1   0   0   0\n")
        f.write("C   0   0   1   0   0\n")
        f.write("U   0   0   0   1   0\n")
        f.write("*   0   0   0   0   1\n")
        f.close()

    def _delete_file(self, filepath):
        """Helper class for deleting a file.

        :type filepath: string
        :param filepath: name (path) of the file
        """
        os.remove(filepath)

    def _generate_plot(self, name):
        """Generates plot using RNAplot. RNAplot must be installed (and present in PATH).

        :type name: string
        :param name: name of the fasta file for which the plot is generated
        """
        # RNAplot cannot plot into different directory -> change dir
        original_dir = os.getcwd()
        os.chdir(self.dir_path)
        # RNAplot always names the file as the fasta header - unfortunately,
        #  fasta header often contains unsupported characters
        # -> create temporary fasta file with header containing only name
        tmp_fasta_path = name + "_tmp.fst"
        tmp_fasta = open(tmp_fasta_path, 'w')
        from_file = open(name + ".fst")
        from_file.readline() # discard first line
        tmp_fasta.write(">" + name + "\n")
        shutil.copyfileobj(from_file, tmp_fasta)
        # for using tmp_fasta as stdin, Python needs to close it & open for reading
        tmp_fasta.close()
        tmp_fasta = open(tmp_fasta_path)
        # perform RNAplot
        args = shlex.split('RNAplot -o svg')
        p = subprocess.Popen(args, stdin=tmp_fasta)
        p.communicate() # waits until the subprocess finishes (and it's safe to remove used file)
        tmp_fasta.close()
        from_file.close()
        self._delete_file(tmp_fasta_path)
        os.chdir(original_dir)

    def _prepare_and_convert_sequences(self):
        """Prepare one big file containing all the secondary structures in a
            ClustalW converted format of all sequences for comparison.
        """
        list_file = open(self.list_filename, 'rU')
        f2 = open(self.all_sequences_filename, 'w')

        for line in list_file:
            line = line.rstrip('\n')
            fst_path = os.path.join(self.dir_path, line) + '.fst'
            if not os.path.isfile(fst_path  ):
                print "Fasta file from the list not found: " + fst_path + "\n";
                print "Fix the file name and run again"
                sys.exit();
            fst_file = open(fst_path)
            f2.write(self._extract_converted_ss(fst_file))

        f2.close()
        list_file.close()

    def _create_output(self, aligned_structures, positions, conservancy_levels):
        """Generates output files - new svg (with "_compare" postfix in name)
            and possibly statistics.

        :type aligned_structures: Dictionary
        :param aligned_structures: list of all aligned structures (with gaps)

        :type positions: List
        :param positions: list of all found positions bordering un/conserved
            regions (basically the (start, stop) Tuples transformed into list)

        :type conservancy_levels: Dictionary
        :param conservancy_levels: list (indexed by position) with conservancy
            levels for each position
        """
        # load expansion segments from file, if supplied
        segments = {}
        if self.segments_file != None:
            with open(os.path.join(self.dir_path, self.segments_file)) as inf:
                segments = eval(inf.read())
        # prepare file for the statistics
        if self.statfile_name != None:
            statFile = open(self.statfile_name, 'a')
            
        conservancyLevelFile = open(os.path.join(self.dir_path, "conservancyLevels" + '_win' + str(self.window_size) + ".txt"), 'w')

        list_file = open(self.list_filename, 'rU')
        for line in list_file:
            name = line.rstrip('\n')
            path = os.path.join(self.dir_path, name)
            conservancyLevelFile.write(name + "\n");

            # generate svg plot
            self._generate_plot(name)

            original_svg = open(path + '_ss.svg')
            colored_svg = open(path + '_win' + str(self.window_size) + '_ss_compare.svg', 'w')

            position = 0
            positionSkipped = 0
            sizeSum = 0
            sizeEsSum = 0
            countEs = 0
            sizeESSum = 0.1
            conservancyLevelString = ""
            for svg_line in original_svg:
                # replace text lines
                if "<text x=\"" in svg_line:
                    # skip gaps (but still count positions in the alignment)
                    while (len(aligned_structures[name]) > position and aligned_structures[name][position] == '-'):
                        position = position + 1
                        positionSkipped = positionSkipped + 1

                    # set default color
                    color = "#333"
                    # if expansion segments are on, determine color by being inside/outside of it
                    if segments.get(name, None) != None:
                        sequence_segments = segments[name]
                        for tt in xrange(len(sequence_segments)):
                            if position - positionSkipped >= sequence_segments[tt][0] and position - positionSkipped <= sequence_segments[tt][1]:
                                color = "blue"
                                if len(conservancy_levels) > position:
                                    sizeESSum = sizeESSum + conservancy_levels[position]
                                break
                    # calculate statistics
                    sizeSum = sizeSum + conservancy_levels[position]
                    # if expansion segments are on & position belongs to expansion segment
                    if color == "blue":
                        sizeEsSum = sizeEsSum + conservancy_levels[position]
                        countEs += 1
                    # rewrite the line with new attributes
                    svg_line = svg_line.replace("<text x=\"", "<text fill=\" " + color + "\" mismatchLvl=\"" + str(conservancy_levels[position]) + "\" pos=\"" + str(position - positionSkipped) + "\" x=\"");
                    if conservancyLevelString == "":
                      conservancyLevelString = str(conservancy_levels[position])
                    else:
                      conservancyLevelString = conservancyLevelString + ", " +  str(conservancy_levels[position])
                    position = position + 1
                elif " height=\"452\" width=\"452\"" in svg_line:
                    svg_line = svg_line.replace(" height=\"452\" width=\"452\"", "")
                colored_svg.write(svg_line)
            # write the statistics
            conservancyLevelString = conservancyLevelString + "\n"
            conservancyLevelFile.write(conservancyLevelString)
            if self.statfile_name != None:
                total_length = position - positionSkipped
                if total_length == 0: # prevent zero division error (in case of empy structure)
                    avgSize = 0
                else:
                    avgSize = sizeSum/total_length
                # complete statistics with expansion segments
                if segments.get(name, None) != None:
                    avgSizeInConserved = sizeEsSum/countEs
                    avgSizeOutConserved = (sizeSum - sizeEsSum)/(total_length - countEs)
                    statFile.write(name + "," + str(self.window_size) + "," + str(sizeSum) + "," + str(sizeESSum) + "," + str(sizeSum/sizeESSum) + "," + str(avgSize) + "," + str(avgSizeInConserved) + "," + str(avgSizeOutConserved) +"\n")
                # simple statistics
                else:
                    statFile.write(name + "," + str(self.window_size) + "," + str(sizeSum) + "," + str(sizeESSum) + "," + str(sizeSum/sizeESSum) + "," + str(avgSize) + ",-,-\n")

        # close everything
        self._delete_file(self.all_sequences_filename)
        list_file.close()
        original_svg.close()
        colored_svg.close()
        conservancyLevelFile.close()
        if self.statfile_name != None:
            statFile.close()

    def compare(self):
        """Perform comparison of all structures named in the list file.
        """
        matrix_path = os.path.join(self.dir_path, "matrix.txt")
        self._create_matrix_file(matrix_path)

        ## PREPARE FILE CONTAINING ALL SEQUENCES IN CONVERTED FORMAT
        self._prepare_and_convert_sequences()

        ## PERFORM MULTIPLE SEQUENCE ALIGNMENT ON ACG-CODED DOT PAREN
        cline = ClustalwCommandline('clustalw2', infile=self.all_sequences_filename, outfile=self.alignment_filename, matrix=matrix_path, type='PROTEIN', gapopen=self.gap_open, gapext=self.gap_extension)
        cline()

        self._delete_file(matrix_path) # matrix is no longer needed

        ## EXTRACT ALL ALIGNED SEQUENCES
        alignment = read_alignment(self.alignment_filename)
        list_file = open(self.list_filename, 'rU')
        aligned_structures = {}

        for line in list_file:
            name = line.rstrip('\n')
            fst_path = os.path.join(self.dir_path, name)
            fst_file = open(fst_path + '.fst')
            header = fst_file.next()
            structure = extract_seq_string_from_aln(alignment, header[1:31].rstrip('\n'))
            # store the structure into a dp file
            # translate the structures back from ClustalW-ACG format to dot-paren
            translated_struct = structure.replace('A', '(').replace('C', ')').replace('G', '.')
            aligned_structures[name] = translated_struct

        list_file.close()

        ## RUN COMPARATION TAGGER ON THE ALIGNMENT
        tagger = ComparationTagger(self.window_size, self.mismatch_percentage, True)

        alignment = read_alignment(self.alignment_filename)
        tagger_result = tagger.find_conserved(alignment)

        conserved_regions = tagger_result[0]
        conservancy_levels = tagger_result[1]

        positions = [] # convert regions starts & stops to a list
        if conserved_regions:
            if isinstance(conserved_regions[0], int):
                print "no conserved regions found"
                return

            for region in conserved_regions:
                positions.append(region[0]);
                positions.append(region[1]);

        self._create_output(aligned_structures, positions, conservancy_levels)
