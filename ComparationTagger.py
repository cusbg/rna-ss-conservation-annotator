#
# ComparationTagger.py
#
# Module for selecting unconserved regions from an alignment.
#
"""This module provides **ComparationTagger** class. A comparation
tagger class is inherited from **ConservationTagger** class
(specifically, from SlidingWindowConservationTagger). The major difference is
that the ComparationTagger considers not only gaps in the alignment, but
actual matches of the characters (which are, for the structure, "(", ")" and ".").
Another difference is that the tagger actually ranks each position with a number.
This number is a percentage of mismatches in a current window, which surrounds
ranked position (e.g. the rank of 100th position for window of size 20
is a percentage of mismatches on positions 80-120).

The first sequence in an alignment is considered a primary one, meaning
other sequences are compared to it. This matter when the first sequence
contains gap, that position is considered mismatch for all other sequences.

The main operation of conservation taggers is through the ``find_conserved``
method. This method ranks all the positions and also creates list of
segments in which the alignment mismatche percentage is above certain limit.

------------------

"""
__author__ = "Jan Pesek"
__contributors__ = "Jan Hajic"
__credits__ = ["Jan Hajic", "David Hoksza", "Josef Panek"]
__maintainer__ = "Jan Pesek"
__email__ = "jpesek89@gmail.com"
__status__ = "Prototype"

from ConservationTagger import SlidingWindowConservationTagger
from collections import Counter

class ComparationTagger(SlidingWindowConservationTagger):
    """Ranks each position of a sequence with conservancy level given by
    amount of mismatches among sequences within a specified window and splits
    sequences into regions divided by certain mismatch percentage in this window.
    """

    def __init__(self, minimum_window=40, gap_percentage=0.3, merge=True):
        """Initializes the ComparationTagger.

        :type minimum_window: int
        :param minimum_window: The number of consecutive positions
            that will be considered for mismatches calculation. Default is 40.

        :type gap_percentage: float
        :param gap_percentage: The maximum percentage of mismatches in a window
            to be considered a conserved region. Default is 0.3.

        :type merge: bool
        :param merge: If set, will merge adjacent and overlapping unconserved
            regions into one.
        """
        super(ComparationTagger, self).__init__()

        self.minimum_window = minimum_window
        self.gap_percentage = gap_percentage
        self.merge = merge

    def _gap_stats(self, alignment, position, mask=None):
        """Computes how many mismatches and matches are in the alignment
            at the given position.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: The alignment to inspect for mismatches.

        :type position: int
        :param position: The position in the alignment at which to count gaps.

        :type mask: list(Boolean)
        :param mask: Not used.

        :rtype: tuple(int, int)
        :returns: A tuple ``(n_mismatches, matches)`` where the members
            are the count of nucleotides and count of mismatches at the given
            position, respectively. The numbers are converted to ``float``
            so that they can be used in mathematical formulae.

        :raises: ValueError
        """

        if (position > len(alignment[0])):
            raise ValueError("Requested match statistics for position > aln. length (position: %d, aln. length: %d)." % (position, len(alignment[0])))

        if position < 0: # starting posiions are calculated from zero position
            position = 0

        column = []
        column = [aln[position] for aln in alignment]

        # find consensus character - the most frequent character on a position
        cnt = Counter()
        for c in column:
            cnt[c] += 1
        consensus = cnt.most_common(1)[0][0]

        # gap cannot be consensuc character
        if consensus == self.ALN_GAP:
            consensus == cnt.most_common(2)[1][0]

        # calculate mismatches
        total_mismatches = len(alignment) - column.count(consensus)
        if len(alignment) == total_mismatches + 1:
            total_mismatches += 1

        return float(len(column) - total_mismatches), float(total_mismatches)

    def find_conserved(self, alignment):
        """Given an alignment, returns a list of ranks of all the alignment
            positions and list of conserved regions.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: A multiple sequence alignment from Bio.Align.

        :rtype: tuple(list(tuple(int,int)), dictionary)
        :returns: A tuple which first element is a list of ``(start, stop)``
             tuples such that all positions between ``start`` and ``stop``
             (including the bounds) are conserved. The regions are sorted
             from lowest to highest position. Note that ``start`` may
             equal ``stop`` if the region is just 1 nucleotide long. Second
             element of the tuple is a dictionary indexed by positions
             containing conservency level (rank) of the position.

             Positions are indexed from 0.
        """
        aln_length = len(alignment[0])
        aln_size = float(len(alignment))

        mismatch_percentage = {}

        window_total = aln_size * (self.minimum_window - 1)

        conserved_regions = []
        # Start building conserved region on window hit. As long
        # as window keeps hitting, extend region by next position.
        current_conserved = [0, 0] # start/stop
        in_conserved_region = False
        for w_start in xrange(aln_length):
            window_end = w_start + (self.minimum_window / 2 - 1)
            # end of the alignment - the window is getting smaller towards the end
            if window_end >= aln_length:
                window_end = aln_length - 1
            non_mismatches, mismatches = self._window_gap_stats(alignment,
                                                                w_start - self.minimum_window / 2,
                                                                window_end
                                                                )
            mismatch_percentage[w_start] = mismatches / window_total
            if mismatch_percentage[w_start] > 1: # can happen in the beginning or end of the structure
                mismatch_percentage[w_start] = 1

            # conserved region indicator - may be used by some other rPredictor components
            if in_conserved_region:
                if mismatches / window_total < self.gap_percentage:
                    continue
                else:
                    current_conserved[1] = w_start + self.minimum_window / 2 - 1
                    # Last window was still good.
                    conserved_regions.append((current_conserved[0],
                                             current_conserved[1]))
                    current_conserved = [0, 0]
                    in_conserved_region = False
            else:
                if mismatches / window_total < self.gap_percentage:
                    current_conserved[0] = w_start - self.minimum_window / 2
                    in_conserved_region = True

        if in_conserved_region:
            current_conserved[1] = aln_length - 1
            conserved_regions.append((current_conserved[0], current_conserved[1]))

        if self.merge:
            conserved_regions = self._merge_regions(conserved_regions)

        return [conserved_regions, mismatch_percentage]
