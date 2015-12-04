#
# ConervationTagger.py
#
# Module for selecting unconserved regions from an alignment.
#
"""This module provides **ConservationTagger** classes. A conservation tagger
is a class that operates on an alignment and finds regions of the alignment
that are (not) conserved well across the aligned sequences. Various methods
of discovering unconserved regions are implemented.

The main operation of conservation taggers is through the ``find_unconserved``
method. This is the method which classes inheriting from
:class:`BaseConservationTagger` have to override in order to do something
useful.

Two conservation tagger algorithms are implemented

* The :class:`SlidingWindowConservationTagger`, which tries to find regions
  such that all its sub-regions of given length contain at least a given
  percentage of gaps

* The :class:`ExtendingWindowConservationTagger`, which tries to find regions
  longer than some minimum such that the overall percentage of gaps in the
  region does not drop below a certain threshold in all subregions sharing
  the maximal unconserved region's start position.

More conservation tagging algorithms planned for the future: most notably,
an HMM model that would take into account more than just gap vs. non-gap
distincions.

------------------

"""
__author__ = "Jan Hajic"
__contributors__ = ""
__credits__ = ["David Hoksza", "Josef Panek", "Michal Klimpera"]
__maintainer__ = "Jan Hajic"
__email__ = "hajicj@ufal.mff.cuni.cz"
__status__ = "Prototype"

import logging


class BaseConservationTagger(object):
    """A ConservationTagger is a class that implements the search for
    conserved and unconserved regions in a Multiple Sequence Alignment.

    Usage:

    >>> alignment = get_alignment('alignment_file.aln')
    >>> ctagger = BaseConservationTagger()
    >>> unconserved_regions = ctagger.find_unconserved(alignment)
    >>> conserved_regions = ctagger.find_conserved(alignment)
    
    etc. The ConservationTagger class is the simplest possible conservation
    tagger: everything is reported as conserved. Subclasses that actually
    do something are

    * ExtendingWindowConservationTagger, which greedily finds windows with
      a given minimum percentage of gaps

    * SlidingWindowTagger, which finds regions such that all their sub-regions
      of given length contain a given minimum percentage of gaps,

    * others, like an HMM tagger (WIP).
    """

    def __init__(self):
        """Initializes the ConservationTagger."""
        self.ALN_GAP = '-'

    def find_unconserved(self, alignment):
        """Given an alignemnt, returns a list of unconserved regions.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: A multiple sequence alignment from Bio.Align.

        :rtype: list(tuple(int,int))
        :returns: A list of ``(start, stop)`` tuples such that all positions
             between ``start`` and ``stop`` (including the bounds) are not
             conserved. The regions are sorted from lowest to highest
             position. Note that ``start`` may equal ``stop`` if the region
             is just 1 nucleotide long.

             Positions are indexed from 0.
        """
        return []

    def find_conserved(self, alignment):
        """Given an alignment, returns a list of conserved regions.
        Complementary to :func:`find_unconserved`; the union of conserved
        and unconserved regions should exactly cover the whole alignment.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: A multiple sequence alignment from Bio.Align.

        :rtype: list(tuple(int,int))
        :returns: A list of ``(start, stop)`` tuples such that all positions
             between ``start`` and ``stop`` (including the bounds) are 
             conserved. The regions are sorted from lowest to highest
             position. Note that ``start`` may equal ``stop`` if the region
             is just 1 nucleotide long.

             Positions are indexed from 0.
        """
        unconserved_regions = self.find_unconserved(alignment)

        if unconserved_regions == []:
            aln_length = len(alignment[0]) 
            return [(0, aln_length-1)]

        conserved_regions = []

        lbound = 0
        rbound = 0
        for u_region in unconserved_regions:
            # First unconserved region starts at overhang
            if u_region[0] == 0:
                lbound = u_region[1] + 1
                continue
            rbound = u_region[0] - 1
            # Unconserved regions are right next to each other
            # (Theoretically possible, although rather nonsensical.
            # Can't be ruled out.)
            if (lbound > rbound):
                lbound = u_region[1] + 1
                continue # There is no conserved region between the
                         # unconserved ones.

            # In the "all right" case:
            conserved_regions.append((lbound, rbound))
            lbound = u_region[1] + 1

        aln_length = len(alignment[0])
        if lbound >= aln_length: # Last unconserved region is overhang
            pass
        else:
            conserved_regions.append((lbound, aln_length-1))

        return conserved_regions

    def _gap_stats(self, alignment, position, mask=None):
        """Computes how many gaps and non-gaps are in the alignment
        at the given position.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: The alignment to inspect for gaps.

        :type position: int
        :param position: The position in the alignment at which to count gaps.

        :type mask: list(Boolean)
        :param mask: If given, will only count the ``i``-th sequence towards
            gap stats if the ``i``-th member of ``mask`` is ``False``.

        :rtype: tuple(int, int)
        :returns: A tuple ``(n_non-gaps, n_gaps)`` where the members
            are the count of nucleotides and count of gaps at the given
            position, respectively. The numbers are converted to ``float``
            so that they can be used in mathematical formulae.

        :raises: ValueError
        """
        
        if (position > len(alignment[0])):
            raise ValueError("Requested gap statistics for position > aln. length (position: %d, aln. length: %d)." % (position, len(alignment[0])))

        column = []
        if mask is None:
            column = [ aln[position] for aln in alignment ]
        else:
            column = [ aln[position] for i, aln in enumerate(alignment) 
                                                        if not mask[i] ]

        total_gaps = 0
        for c in column:
            if c == self.ALN_GAP:
                total_gaps += 1

        return float(len(column) - total_gaps), float(total_gaps)

    def _window_gap_stats(self, alignment, start, stop, mask=None):
        """Computes how many gaps and non-gaps are in the alignment
        in the ``(start, stop)`` region (incl. both ``start`` and ``stop``
        positions).

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: The alignment to inspect for gaps.

        :type start: int
        :param start: The first position at which to count gap statistics

        :type stop: int
        :param stop: The last position at which to count gap statistics.

        :type mask: list(list(Boolean))
        :param mask: If given, will only count at each position in the window
            the ``i``-th sequence towards gap stats if the ``i``-th member of 
            ``mask`` is ``False``.

        :rtype: tuple(int, int)
        :returns: A tuple where the first member is the number of non-gap
            characters in the given region of the alignment and the second is
            the number of gap characters.

        :raises: ValueError
        """
        # Sanity check
        if (stop < start):
            raise ValueError("Requested gap statistics for window with stop < start (start: %d, stop: %d)." % (start, stop))
        if (stop > len(alignment[0])):
            raise ValueError("Requested gap statistics for window with stop > aln. length (stop: %d, aln. length: %d)." % (stop, len(alignment[0])))

        non_gaps = 0.
        gaps = 0.
        for i, position in enumerate(xrange(start, stop+1)):

            position_mask = None
            if mask is not None:
                position_mask = [ m[i] for m in mask ]
           
            p_non_gaps, p_gaps = self._gap_stats(alignment, position,
                                                 position_mask)
            non_gaps += p_non_gaps
            gaps += p_gaps

        return (non_gaps, gaps)

    def _shave_edges(self, alignment, region):
        """Removes from the region's edges all columns with no gaps.

        :type alignment: Bio.Align.Alignment
        :param alignment: The alignment with respect to which we are performing
            the operation.

        :type region: tuple(int, int)
        :param region: A (start, stop) delimitation of the region, off of which
            we want to "shave" positions with no gaps.

        :rtype: tuple(int, int)
        :returns: A new (start, stop) region that is the largest sub-region
            of the ``region`` parameter such that there are gaps in the alignment
            at both the ``start`` and ``stop`` position. (May be equivalent to the
            input region.)
        """
        left = region[0]
        right = region[1]

        l_non_gaps, l_gaps = self._gap_stats(alignment, left)
        while l_gaps == 0:
            left += 1
            l_non_gaps, l_gaps = self._gap_stats(alignment, left)

        r_non_gaps, r_gaps = self._gap_stats(alignment, right)
        while r_gaps == 0:
            right -= 1
            r_non_gaps, r_gaps = self._gap_stats(alignment, right)

        return (left, right)


class SlidingWindowConservationTagger(BaseConservationTagger):
    """Searches for continuous sequences of windows with a given minimum
    percentage of alignment gaps."""

    def __init__(self, minimum_window=1, gap_percentage=0.3, merge=True,
                 shave=True, disregard_ends=True):
        """Initializes the ConservationTagger.

        :type minimum_window: int
        :param minimum_window: The minimum number of consecutive positions
            that will be considered for an unconserved region. Default is 1.

        :type gap_percentage: float
        :param gap_percentage: The minimum percentage of gaps in a window
            to be considered an unconserved region. Default is 0.3. You will
            want to estimate this number from your data at some point later
            on.

        :type merge: bool
        :param merge: If set, will merge adjacent and overlapping unconserved
            regions into one.

        :type shave: bool
        :param shave: If set, will shave off each unconserved region all
            columns at the start or end of the region which don't have a gap.

        :type disregard_ends: bool
        :param disregard_ends: If given, will not count towards the gap
            threshold gaps that are at the starts and ends of sequences.
        """
        super(SlidingWindowConservationTagger, self).__init__()

        self.minimum_window = minimum_window
        self.gap_percentage = gap_percentage
        self.merge = merge
        self.shave = shave
        self.disregard_ends = disregard_ends

    def find_unconserved(self, alignment):
        """Given an alignemnt, returns a list of unconserved regions.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: A multiple sequence alignment from Bio.Align.

        :rtype: list(tuple(int,int))
        :returns: A list of ``(start, stop)`` tuples such that all positions
             between ``start`` and ``stop`` (including the bounds) are not
             conserved. The regions are sorted from lowest to highest
             position. Note that ``start`` may equal ``stop`` if the region
             is just 1 nucleotide long.

             Positions are indexed from 0.
        """
        aln_length = len(alignment[0])
        aln_size = float(len(alignment))

        window_total = aln_size * self.minimum_window

        unconserved_regions = []

        alignment_ends_mask = None
        if self.disregard_ends:
            alignment_ends_mask = self._get_aln_ends_mask(alignment)

        # Start building unconserved region on window hit. As long
        # as window keeps hitting, extend region by next position.
        current_unconserved = [0,0] # start/stop
        in_unconserved_region = False
        for w_start in xrange(aln_length - self.minimum_window + 1):

            if w_start % 100 == 0:
                logging.debug('  Tagging at position %d' % w_start)

            # Implemented for the disregard_ends option.
            window_mask = None
            if alignment_ends_mask is not None:
                window_mask = [ a[w_start:w_start + self.minimum_window] 
                                for a in alignment_ends_mask ]

            non_gaps, gaps = self._window_gap_stats(alignment, w_start,
                                       w_start + (self.minimum_window - 1),
                                       window_mask)
            
            if in_unconserved_region:
                if gaps / window_total > self.gap_percentage:
                    continue
                else:
                    current_unconserved[1] = w_start + self.minimum_window - 1
                      # Last window was still good.
                    unconserved_regions.append((current_unconserved[0], 
                                                current_unconserved[1]))
                    current_unconserved = [0,0]
                    in_unconserved_region = False
            else:
                if gaps / window_total > self.gap_percentage:
                    current_unconserved[0] = w_start
                    in_unconserved_region = True

        if in_unconserved_region:
            current_unconserved[1] = aln_length - 1

            unconserved_regions.append((current_unconserved[0],
                                        current_unconserved[1]))

        if self.merge:
            unconserved_regions = self._merge_regions(unconserved_regions)

        if self.shave:
            unconserved_regions = [ self._shave_edges(alignment, region) 
                                    for region in unconserved_regions ]

        return unconserved_regions

    def _merge_regions(self, regions):
        """Merges all overlapping and adjacent regions into one.

        :type regions: list(tuple(int,int))
        :param regions: A list of regions that should be checked for overlap
            and adjacency.

        :rtype: list(tuple(int,int))
        :returns: A list of the merged regions.
        
        """
        if len(regions) == 0:
            return []
        output_regions = []
        current_output = list(regions[0])
        for region in regions[1:]:
            if current_output[1] + 1 >= region[0]: # Adjacent or overlap
                current_output[1] = region[1]
            else:
                output_regions.append(tuple(current_output))
                current_output = list(region)
        output_regions.append(tuple(current_output))
        return output_regions

    def _get_aln_ends_mask(self, alignment):
        """Returns a mask such that it is True for gaps at ends of sequences
        and False everywhere else.

        >>> CT._get_aln_ends_mask(['--AA--', '---UUA', 'G-C-A-'])
        [[1, 1, 0, 0, 1, 1], [1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1]]

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: The multiple sequence alignment on which the
            conservation tagger operates.

        :rtype: list(list(Boolean))
        :returns: A mask for the entire alignment such that True entries are
            for sequence regions at starts and ends before the first non-gap
            is found.
        """
        aln_ends_mask = [ [] for a in alignment ]

        for i, sequence in enumerate(alignment):

            current_mask = [ False for s in sequence ]

            for position in xrange(len(sequence)):
                if sequence[position] != self.ALN_GAP:
                    break
                else:
                    current_mask[position] = True

            for position in xrange(len(sequence) - 1, 0, -1):
                if sequence[position] != self.ALN_GAP:
                    break
                else:
                    current_mask[position] = True

            aln_ends_mask[i] = current_mask

        return aln_ends_mask


class ExtendingWindowConservationTagger(BaseConservationTagger):
    """Searches for maximal windows with a given minimum percentage
    of alignment gaps."""

    def __init__(self, minimum_window=1, gap_percentage=0.3, 
                 cutback=True):
        """Initializes the ConservationTagger.

        :type minimum_window: int
        :param minimum_window: The minimum number of consecutive positions
            that will be considered for an unconserved region. Default is 1.

        :type gap_percentage: float
        :param gap_percentage: The minimum percentage of gaps in a window
            to be considered an unconserved region. Default is 0.3. You will
            want to estimate this number from your data at some point later
            on.

        :type cutback: bool
        :param cutback: If set, the tagger will go backward from each
            unconserved region's end, find an unconserved window in much the
            same way as the forward pass and the output region is limited
            by the end of the backward window.
        """
        super(ExtendingWindowConservationTagger, self).__init__()

        self.minimum_window = minimum_window
        self.gap_percentage = gap_percentage
        self.cutback = cutback

    def find_unconserved(self, alignment):
        """Given an alignemnt, returns a list of unconserved regions.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: A multiple sequence alignment from Bio.Align.

        :rtype: list(tuple(int,int))
        :returns: A list of ``(start, stop)`` tuples such that all positions
             between ``start`` and ``stop`` (including the bounds) are not
             conserved. The regions are sorted from lowest to highest
             position. Note that ``start`` may equal ``stop`` if the region
             is just 1 nucleotide long.

             Positions are indexed from 0.
        """
        aln_length = len(alignment[0])
        aln_size = float(len(alignment))
          # We'll be dividing things by this number.

        unconserved_regions = []

        w_start = 0
        while w_start < aln_length - self.minimum_window + 1:
            non_gaps, gaps = self._window_gap_stats(alignment, w_start, 
                                          w_start + (self.minimum_window - 1))

            logging.debug('At w_start %d, gaps %d, non-gaps %d' % (w_start, gaps, non_gaps))

            if gaps / (gaps + non_gaps) >= self.gap_percentage:
                test_pos = w_start + self.minimum_window
                while ((gaps / (gaps + non_gaps) >= self.gap_percentage) 
                       and (test_pos < aln_length)):
                    p_non_gaps, p_gaps = self._gap_stats(alignment, test_pos)
                    non_gaps += p_non_gaps
                    gaps += p_gaps
                    test_pos += 1
                # When the cycle ends, either the percentage dropped below
                # minimum when taking test_pos into account, or test_pos is
                # at aln_length.
                w_end = test_pos - 1
                if self.cutback:
                    # Perform backward pass
                    # Find intersection
                    w_end = self._cutback_pass(alignment, w_end)
                
                unconserved_regions.append((w_start, w_end))

                # After finding an unconserved region, start looking only
                # to the right of it.
                w_start = w_end + 1
            else:
                w_start += 1

        return unconserved_regions

    def _cutback_pass(self, alignment, w_end):
        """Runs the cutback pass of Extending Window.
        We only want to discover a single window. Returns """
        while w_end > self.minimum_window - 1:
            non_gaps, gaps = self._window_gap_stats(alignment, 
                                     w_end - self.minimum_window + 1, w_end)
            if gaps / (non_gaps + gaps) >= self.gap_percentage:
                return w_end 
            else:
                w_end -= 1

        raise ValueError('For some reason, _cutback_pass hasn\'t found a backward window.')


class OptimalTagger(BaseConservationTagger):
    """This tagger finds the optimal conservation assignment given
    an alignment and a structure of one of the molecules in the alignment.

    Optimality is defined as the assignment cf conserved regions that
    maximizes the f-score between gap recall in unconserved regions and
    base pair recall in conserved regions. Gap recall is defined as::

      # of gaps in unconserved / total # of gaps

    and base pair recall is defined as::

      # of BPs with both ends in conserved / total # of BPs
    """
    def __init__(self):
        """Initializes the tagger."""
        super(self, OptimalTagger).__init__()

    def find_unconserved(self, alignment, template, target):
        """Given an alignment, a template structure and a target structure,
         finds the optimal set of unconserved regions so that prediction
         precision and recall is maximized. This tagger finds an upper bound
         for a given template-target pair.

        :type alignment: Bio.Align.MultipleSeqAlignment
        :param alignment: A multiple sequence alignment from Bio.Align.

        :type template: Secstruc.PseudoknotSecstruc
        :param template: The template structure for prediction.

        :type target: Secstruc.PseudoknotSecstruc
        :param target: The target structure for prediction.

        :rtype: list(tuple(int,int))
        :returns: A list of ``(start, stop)`` tuples such that all positions
            between ``start`` and ``stop`` (including the bounds) are not
            conserved. The regions are sorted from lowest to highest
            position. Note that ``start`` may equal ``stop`` if the region
            is just 1 nucleotide long.

            Positions are indexed from 0.

        """
        # How to find the optimal set of unconserved regions?
        #  - find areas where most errors are

        #
        # If we do not allow outpairs:
        #  - make an inclusion structure, mark # of base pairs (top-down)
        #    we lose by accepting each BP span
        #  - mark # of gaps in each BP span (bottom-up)
        #  - What about unconserved loop sections in junctions?
        #
        return []
