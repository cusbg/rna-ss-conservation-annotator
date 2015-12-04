#!/usr/bin/env python
#
# Secstruc.py
"""Module for handling RNA secondary structures built around dot-bracket
notation.

The clases provide methods for discovering *structural features*, certain
well-defined substructures that serve as building blocks of the secondary
structure.

------------------------------

"""
__author__ = "Kristian Rother"
__contributors__ = "Tomasz Puton, Jan Hajic"
__credits__ = ["Lukasz Kozlowski, "
               "Natalia Szostak, "
               "Joanna Kasprzak, "
               "Sandra Smit"]
__maintainer__ = "Jan Hajic"
__email__ = "hajicj@ufal.mff.cuni.cz"
__status__ = "Development"

import logging
import re

import PKnotSecstrucConstants as Constants
from PKnotResolver import PKnotResolver


class SecstrucError(Exception):
    """Raised in case of secondary structure error not covered by some
    more specific exception."""
    pass


class PseudoknotError(Exception):
    """Raised if a secondary structure encounters unexpected pseudoknots."""
    pass


class Secstruc(object):
    """Secondary structure string with indices for each position.

    Can handle pieces of secondary structure which are not continuous.
    However, when indexing, all ``__getitem__``-based operations
    are using 0-based indexing, both into the structural dot-parent string
    and into the array of "real" indexes.

    Note that secondary structures with base pairs leading outside the structure
    are allowed, as sometimes molecules form base pairs outside themselves. To
    check that the structure does *not* contain such "half-pairs", use the
    :meth:`validate` method.
    """
    def __init__(self, secstr, indices=None, brackets=('(', ')'),
                 validate=False):
        """Initializes a secondary structure object.

        :type secstr: string
        :param secstr: A dot-paren secondary structure string.

        :type indices: list(int)
        :param indices: A list of indices that corresponds to positions in the
            ``secstr`` parameter. These indices are NOT used for accessing the
            secondary structure using the ``__getitem__`` method - they do not
            necessarily have to start with 0, or even be ordered (but that would
            be a *bad* idea).

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :type validate: bool
        :param validate: If set, will raise ``SecstrucError`` if the given
            structure contains an unopened/unclosed base pair.
        """
        self.secstruc = secstr
        if indices is None:
            indices = range(len(self))
        if len(self.secstruc) != len(indices):
            raise SecstrucError("Cannot create Secstruc object (%s %s):" % (
                secstr, str(indices)) +
                                " db string and indices lengths don't match." )
        self.indices = indices
        self.indices = indices

        self.base_pairs = self.get_base_pair_indices_list(brackets)
        self._base_pair_index = set(self.base_pairs)

        if validate:
            if not self.validate():
                raise SecstrucError("Cannot create Secstruc object "
                                    "(%s %s): " % (secstr, str(indices)) +
                                    "validation failed.")

    def find_overhang5(self, brackets=('(', ')')):
        """Returns overhang on the 5' end as Secstruc objects.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(Secstruc)
        :returns: The secondary structure object corresponding to the
            5'-overhang sliced from the structure wrapped into a list (although
            there is only one overhang of each type, of course; this is to make
            the processing of structural elements easier).
        """
        b0 = brackets[0]

        n = len(self)
        # "No overhang" case (guards against None return)
        if self.secstruc[0] != '.':
            return [self[:1]]

        overhang = []
        for i in xrange(2, n):
            char = self.secstruc[i]
            if char == '.':
                continue
            elif char == b0:
                overhang = [self[:i+1]]
                break
            else:
                logging.warn('Found unexpected character when searching for'
                             '5-overhang: %s at position %d,' % (char, i) +
                             ' empty overhang returned.')
                break

        logging.debug('Returning 5-overhang: %s' % str(overhang))
        return overhang

    def find_overhang3(self, brackets=('(', ')')):
        """Returns overhang on the 3' end as Secstruc objects.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(Secstruc)
        :returns: The secondary structure object corresponding to the
            3'-overhang sliced from the structure wrapped into a list (although
            there is only one overhang of each type, of course; this is to make
            the processing of structural elements easier).
        """
        b1 = brackets[1]

        n = len(self)

        # "No overhang" case (guards against None return)
        if self.secstruc[n-1] != '.':
            logging.debug('Last character in string not an unpaired base:'
                          ' %s' % self.secstruc[n-1])
            logging.debug('self.secstruc = %s' % self.secstruc)
            return [self[n-1:]]

        overhang = []
        for i in xrange(n-1, 0, -1):
            char = self.secstruc[i]
            if char == '.':
                continue
            elif char == b1:
                overhang = [self[i:]]
                logging.debug('Found 3\' overhang end %s at %d.' % (char, i))
                break
            else:
                logging.warn('Found unexpected character when searching for'
                             ' 3-overhang: %s at position %d,' % (char, i) +
                             ' empty overhang returned.' )
                break

        logging.debug('Returning 3-overhang: %s' % str(overhang))
        return overhang

    def find_bulges(self, brackets=('(', ')')):
        """Returns a list of all bulge structural elements.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(Secstruc)
        :returns: A list of the secondary structure object corresponding to the
            bulges sliced from the structure.
        """
        bulges = []
        for i, j in self.base_pairs:
            bulge, ind = self._check_bulge(i, j, brackets)
            if bulge:
                bulges.append(bulge)
        return bulges

    def find_helices(self):
        """Returns a list of all helices in the structure.

        A helix is defined as a set of two regions where the i-th residue
        from the 5'-end of one is paired to the i-th residue from the 3'-end
        of the other. This is a helix::

          ((()))

        This is not a helix, but two helices::

          ((.()))

        :rtype: list(Secstruc)
        :returns: A list of Secstrucs corresponding to helices in the structure.

        """
        # return self.extract_elements('helix')
        helices = []
        begin = 0
        end = 0
        length = 0

        for i, j in self.get_base_pair_indices():
            # check if old helix continues
            if i == begin + length and j == end - length:
                length += 1
            else:
                # save old helix
                if length > 0:
                    helices.append(self[begin:begin+length] +
                                   self[end-length+1:end+1])
                # new helix starts
                begin = i
                end = j
                length = 1
        if length > 0:
            helices.append(self[begin:begin+length] + self[end-length+1:end+1])

        return helices

    def find_hairpins(self, brackets=('(', ')')):
        """Returns a list of all hairpin structural elements.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(Secstruc)
        :returns: A list of the secondary structure object corresponding to the
            bulges sliced from the structure.
        """
        loops = []
        for i, j in self.base_pairs:
            loop = self._check_hairpin(i, j, brackets)
            if loop:
                loops.append(loop)
        return loops

    def find_internal_loops(self, brackets=('(', ')')):
        """Returns a list of all internal loop structural elements.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(Secstruc)
        :returns: A list of the secondary structure object corresponding to the
            internal loops sliced from the structure.
        """
        loops = []
        for i,j in self.base_pairs:
            loop = self._check_internal_loop(i, j, brackets)
            if loop:
                loops.append(loop)
        return loops

    def find_junctions(self, brackets=('(', ')')):
        """Returns a list of all junction structural elements.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(Secstruc)
        :returns: A list of the secondary structure object corresponding to the
            junctions sliced from the structure.
        """
        b0 = brackets[0]
        b1 = brackets[1]

        junctions = []
        for i, j in self.get_base_pair_indices():
            # check if i,j enclose a junction.
            nested_bps = 0
            skip_helix = 0
            junction = self[i]
            k = i+1
            while k < j:
                if self.secstruc[k] == b1:
                    skip_helix -= 1
                    if skip_helix == 0:
                        nested_bps += 1
                if skip_helix == 0:
                    junction += self[k]
                if self.secstruc[k] == b0:
                    skip_helix += 1
                k += 1
            junction += self[j]
            if nested_bps >= 2:
                junctions.append(junction)
        return junctions

    def _check_hairpin(self, i, j, brackets=('(', ')')):
        """Checks whether the given region ``(i, j)`` is a hairpin loop.

        .. warning::

          Currently undefined for hairpins with uncontignuous indices.

        :type i: int
        :param i: The index of the first residue to check. Index into the
            secstruc string, NOT an index kept in the ``indices`` array!

        :type j: int
        :param j: The index of the last residue to check. Index into the
            secstruc string, NOT an index kept in the ``indices`` array!

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: Boolean
        :returns: True if given region of the structure is a hairpin loop,
            False if not.
        """
        idx = i + 1
        # See check_bulge for how to adapt for dealing with PKs.
        while self.secstruc[idx] not in brackets:  # j is ), so idx < j
                                                   # not needed.
            idx += 1
        if idx == j: # Found loop!
            return self[i:j+1]  # Slice fencepost...
        else:
            return None

    def _check_internal_loop(self, i, j, brackets=('(', ')')):
        """Checks whether the given region ``(i, j)`` has an internal loop.
        This is NOT the same as checking for a hairpin, where we check that the
        entire region is a hairpin loop. Here, we check whether a base pair
        ``(i,j)`` is the outer base pair of an internal loop::

            i                             j
            (....( /* something */ )......)

        We look for an inner base pair.

        :type i: int
        :param i: The index of the 5'-end of the outer base pair.

        :type j: int
        :param j: The index of the 3'-end of the outer base pair.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: Boolean
        :returns: True if given region of the structure has an internal loop,
            False if not.        
        """
        
        b0 = brackets[0]
        b1 = brackets[1]

        idx_i = i+1
        while (self.secstruc[idx_i] == '.') and (idx_i < j-1):
            idx_i += 1
        if (self.secstruc[idx_i] == b0) and (idx_i - i > 1):
            idx_j = j-1
            while (self.secstruc[idx_j] == '.') and (idx_j > idx_i):
                idx_j -= 1
            if (self.secstruc[idx_j] == b1) and (j - idx_j > 1) and \
               ((idx_i, idx_j) in self._base_pair_index):  # against (..()()...)
                return self[i:idx_i+1] + self[idx_j:j+1]
        return None

    def _check_bulge(self, i, j, brackets=('(', ')')):
        """Checks whether the given region ``(i, j)`` has a bulge.
        This is NOT the same as checking for a hairpin, where we check that the
        entire region is a hairpin loop. Here, we check whether a base pair
        ``(i,j)`` is the outer base pair of a bulge::

            i                     j
            (..( /* something */ ))

        We look for an inner base pair that is adjacent at one end and
        non-adjacent at the other.

        :type i: int
        :param i: The index of the 5'-end of the outer base pair.

        :type j: int
        :param j: The index of the 3'-end of the outer base pair.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: Boolean
        :returns: True if given region of the structure has a bulge, False if
            not.        
        """

        if i > len(self) - 3 or j == 0:
            return None, None
        if j - i < 3:
            return None, None

        b0 = brackets[0]
        b1 = brackets[1]

        # In order to incorporate pseudoknots, change checking for '.' to
        # checking for '.', ']' or '[' when extending bulged side.

        # In order to exclude pseudoknotted items from the bulge,
        # do not append/prepend given position to the bulged side if
        # the current nucleotide is pseudoknotted.

        # Only need to check for bulges is the structure is (. )) or (( .)
        if (self.secstruc[i+1] not in brackets) and (self.secstruc[j-1] == b1):
            # Checks the i-side bulge. 
            j_side = self[j-1:j+1] # Slice notation fencepost peculiarity.
            i_side = self[i]
            i_side_idx = i+1
            # Push all unpaired bases adjoining position i on the i-strand
            # to the bulge.
            # Currently implemented variant: if it hits a pseudoknot, doesn't
            # report a bulge.
            while (i_side_idx + 1 < j - 1) and (self.secstruc[i_side_idx] not in brackets):
                i_side += self[i_side_idx]
                i_side_idx += 1
            # If we hit j and there were all dots: this wasn't a base pair! 
            # We may have a bulge! Need to check for junction: (..()())
            if self.secstruc[i_side_idx] == b0:
                if (i_side_idx, j-1) not in self._base_pair_index:
                    return None, None # False alarm - junction.
                i_side += self[i_side_idx]
                bulge = i_side + j_side
                ### DEBUG
                #print 'Found bulge. i_side=' + i_side.secstruc + ', j_side=' + j_side.secstruc
                return bulge, i+1
            else:
                return None, None

        # If we don't have a bulge: check the j-side bulge
        elif (self.secstruc[i+1] == b0) and (self.secstruc[j-1] not in brackets):
            # Checks the j-side bulge.
            i_side = self[i:i+2] # Slice notation fencepost peculiarity.
            j_side = self[j]
            j_side_idx = j - 1
            while (j_side_idx > i + 1) \
                    and (self.secstruc[j_side_idx] not in brackets):
                j_side = self[j_side_idx] + j_side  # Prepend!
                j_side_idx -= 1

            if self.secstruc[j_side_idx] == b1:
                # We may have a bulge! Need to check for junction: (()()..)
                if (i+1, j_side_idx) not in self._base_pair_index:
                    return None, None  # False alarm - junction.
                j_side = self[j_side_idx] + j_side
                bulge = i_side + j_side
                return bulge, j-1
            else:
                return None, None
        else:  # Whatever else happens, NOT A BULGE.
            return None,None

    def get_base_pair_indices(self, brackets=('(', ')')):
        """Generates (i,j) base pair tuples from the dot-paren string given in
        the init structure.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: Generator of tuples
        :returns: Yields a list of the base pairs represented as (start, end)
            indices.
        """
        n = len(self)
        i = 0
        b0 = brackets[0]
        b1 = brackets[1]
        while i < n - 1:
            # find next open pair
            while i < n - 1 and self.secstruc[i] != b0: i += 1
            # find corresponding pair
            j = i + 1
            bps = 1
            while j < n and bps > 0:
                if self.secstruc[j] == b0:
                    bps += 1
                if self.secstruc[j] == b1:
                    bps -= 1
                if bps == 0:
                    yield i, j
                j += 1
            i += 1

    def get_base_pair_indices_list(self, brackets=('(', ')')):
        """Generates (i,j) base pair tuples from the dot-paren string given in
        the init structure.

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs.

        :rtype: list(tuple(int,int))
        :returns: Returns a list of the base pairs represented as (start, end)
            indices.
        """
        bps = []
        for i, j in self.get_base_pair_indices(brackets):
            bps.append((i, j))
        return bps

    def validate(self, diagnostics=False, brackets=('(', ')')):
        """Checks whether all base pairs opened in the structure are also
        closed.

        :type diagnostics: bool
        :param diagnostics: If set, if validation fails, will additionally
            check for exact reason why validation failed and find the offending
            positions in the structure. (Currently not implemented.)

        :type brackets: tuple(char, char)
        :param brackets: The brackets used in the ``secstr`` dot-paren string
            to denote base pairs. (The function only checks for one level
            of psueoknottedness, as identified by bracket type.)

        :rtype: bool
        :returns: True if structure is valid (no base pairs left unclosed or
            unopened), False otherwise.
        """
        n_bps = len(self.base_pairs)
        n_open = 0
        n_closed = 0
        for c in self.secstruc:
            if c == brackets[0]:
                n_open += 1
            if c == brackets[1]:
                n_closed += 1

        if n_open == n_closed:
            if n_bps != n_open:
                logging.debug('All pairs at bracket level {0} match, but there'
                              ' are other levels.'.format(brackets[0]))
            return True

        else:
            logging.debug('Invalid secstruc: level {0}, n_bps={1}, n_open={2},'
                          ' n_closed={3}'.format(brackets[0],
                                                 n_bps, n_open, n_closed))

            if diagnostics:
                self._validation_diagnostics()

            return False

    def _validation_diagnostics(self):
        """Tries to diagnose what kind of problem is there with an invalid
        structure.

        NOT IMPLEMENTED."""
        raise NotImplementedError()

    def filter_outpairs(self):
        """Creates a structure from itself that has unpaired nucleotides
        instead of outpair ends.

        >>> structure = PseudoknotSecstruc('..)((..[[)[..)]](]..((')
        >>> structure.filter_outpairs().secstruc
        '...((..[[)[..)]].]....'

        :returns: A structure of the same class with the outpair ends removed
            from dot-paren representation.
        """
        new_dp = make_dotbracket_from_bplist(len(self), self.base_pairs)
        original_dp_chars = [ '.' for _ in new_dp ]
        for i, char in enumerate(new_dp):
            if char == '.':
                original_dp_chars[i] = '.'
            else:
                original_dp_chars[i] = self.secstruc[i]

        original_dp_string = ''.join(original_dp_chars)
        filtered_struct = self.__class__(original_dp_string, self.indices)

        # The filtered structure may have different results for pseudoknot
        # levels. We need to revert that to the original structure.
        return filtered_struct

    def find_elements(self, element_type):
        """Finds elements of a given type. The permitted types:

        * 3'-overhang
        
        * 5'-overhang

        * bulge

        * hairpin

        * helix

        * internal loop

        * junction

        .. warning::

          This method is subject to deprecation.

        :type element_type: str
        :param element_type: The element type to search for.

        :rtype: list(Secstruc)
        :returns: A list of the secondary structure object corresponding to the
            structural elements given by ``element_type`` sliced from the
            structure.        
        """
        result = []
        n = len(self)
        i = 0
        while i < n-1:
            j = n
            while j > i:
                ele = self[i:j]
                if ele.get_type() == element_type:
                    result.append(ele)
                j -= 1
            i += 1
        return result

    def find_Secstruc_elements(self):
        """
        Takes a secondary structure, and returns a list of
        component secondary structure elements.

        .. warning::

            Deprecated.

        """

        logging.warn('This method is deprecated. Please refer to methods '
                     'find_helices(), find_hairpins(), etc.')

        result = []
        overhang = self.find_overhang5()
        if overhang: result.append(overhang)
        overhang = self.find_overhang3()
        if overhang: result.append(overhang)
        loops = self.find_loops()
        helices = self.find_helices()
        bulges = self.find_bulges()
        junctions = self.find_junctions()
        result += loops + bulges + helices + junctions
        return result

    def get_type(self):
        """Returns a type classification of this secondary structure as
        a string, if it is some basic type.

        .. warning::

            Deprecated.

        """
        if re.search('^\.+\($',self.secstruc): return "5'-overhang"
        if re.search('^\)\.+$',self.secstruc): return "3'-overhang"
        if re.search('^\(\(+\)\)+$',self.secstruc):
            # check both stems have equal length
            half = len(self.secstruc)/2
            if len(self.secstruc)%2 == 0 \
               and re.search('^\(+$',self.secstruc[:half]) \
               and re.search('^\)+$',self.secstruc[half:]):
                return 'helix'
            else: return None
        if re.search('^\(\.\.+\)$',self.secstruc): return "loop"
        if re.search('^\(\(\)\.+\)|\(\.+\(\)\)$',self.secstruc): return "bulge"
        if re.search('^\((.*\(\)\.*)+\)$',self.secstruc): return 'junction'
        return None

    def __eq__(self, other):
        """Compares two secondary structure objects. Equality is defined if both
        the secondary structure dot-paren strings and the indices array are the
        same.

        :type other: Secstruc
        :param other: The secondary structure to which this one is compared.

        :rtype: Boolean
        :returns: True if the structures are the same, false otherwise.
        """

        # DEBUG
        #print "Comparing: ", self, "is", type(self), "//", other, "is", type(other)

        if not isinstance(other, Secstruc):
            other = self.__class__(other)
            # DEBUG
            #print "  converting", other, "to secstruc(other)"
            #print "  self.secstruc =", self.Secstruc, "// other.Secstruc =", other.Secstruc
            #print "  self.indices =", self.indices, "// other.indices =", other.indices
        if self.secstruc == other.secstruc \
           and self.indices == other.indices:
            return True
        
    def __add__(self, other):
        """Concatenates the secondary structures together.

        :type other: Secstruc
        :param other: The secondary structure to concatenate to this one.

        :rtype: Secstruc
        :returns: A new secondary structure obtained as a concatenation of
            the original dot-paren strings and indices. Note that this may mean
            duplicate indices, if care is not taken to avoid the situation! This
            method is intended for concatenating substructures obtained by the
            ``__getitem__`` method from **one structure** into a noncontignuous
            substructure (like a helix or junction structural element).

            .. note::

                Duplicate indices do *not* obstruct :meth:`__getitem__` action,
                because it does not index the structure using the ``indices``
                attribute.
        """
        newsec = self.secstruc + other.secstruc
        newind = self.indices + other.indices
        output = self.__class__(newsec, newind)
        #logging.debug('Result of adding\n%s\nand\n%s\nis\n%s' % (self,
        #                                                         other,
        #                                                         output))
        return output

    def __getitem__(self, index):
        """Returns a substructure: either a slice or a single position. 

        :type index: int or slice
        :param index: The index or indices of the substructure. A request for
            the ``i``-th item returns the ``i``-th member of the structure's
            ``secstruc`` string and the ``indices`` array, regardless of whether
            ``indices[i] == i``.            

        :rtype: Secstruc
        :returns: A substructure - either of length 1, if a single index is
            given, or a region, if a slice is given.
        """
        if type(index) == slice:
            newsec = self.secstruc[index.start:index.stop]
            newind = self.indices[index.start:index.stop]
            output =  self.__class__(newsec, newind)
            return output
        else:
            output = self.__class__(self.secstruc[index], [self.indices[index]])
            return output

    def __repr__(self):
        return self.secstruc + ';' + str(self.indices)

    def __len__(self):
        return len(self.secstruc)

    def __str__(self):
        return self.secstruc


class PseudoknotSecstruc(Secstruc):
    """A class for secondary structures that contain pseudoknots.

    Uses constants from :class:`PKnotSecstrucConstants` to abstract from
    the round brackets used in :class:`Secstruc` to denote base pairs.
    """

    def get_base_pair_indices(self, brackets=None):
        """Generates (i,j) base pair tuples from the dot-paren string given in
        the init structure.

        .. note::

          The algorithm used here is not as efficient as it could be (searches
          forward for the closing symbol instead of keeping stacks of open
          symbols).

        :type brackets: tuple(char, char)
        :param brackets: Ignored in this method of PseudoknotSecstruc - a
            different algorithm is used for retrieving base pairs.

        :rtype: Generator of tuples
        :returns: Yields a list of the base pairs represented as (start, end)
            indices.
        """

        n = len(self)
        i = 0

        # DEBUG
        #logging.debug('Running pseudoknotted get_base_pair_indices...')
        #print "Open symbols:", PseudoknotSecstruc.open_symbols
        #print "len(self) =", n

        while i < n - 1:
            # find next open pair
            while ((i < n - 1) and (str(self.secstruc[i]) not in
                                  Constants.open_symbols)):
                # DEBUG
                #print "Searching for open_symbol, at", i, "| symbol", \
                # self.secstruc[i], "| open_symbols =", \
                # PseudoknotSecstruc.open_symbols
                #print (i<n-1 and (str(self.secstruc[i]) not in \
                #  PseudoknotSecstruc.open_symbols))
                #print type(PseudoknotSecstruc.open_symbols[0]),
                #           type(self.secstruc[i])

                i += 1

            # Check against 3'-end
            if not (i < n - 1):
                break

            # find corresponding pair
            open_symbol = str(self.secstruc[i])

            # DEBUG
            #print "self.secstruc[i] is", type(self.secstruc[i]), \
            #      "self.secstruc is", type(self.secstruc)
            #print "open_symbol =", open_symbol, "/ is", type(open_symbol)

            close_symbol = Constants.symbols[open_symbol]

            j = i + 1
            bps = 1
            while j < n and bps > 0:
                if str(self.secstruc[j]) == open_symbol:
                    bps += 1
                if str(self.secstruc[j]) == close_symbol:
                    bps -= 1
                if bps == 0:
                    yield i, j
                j += 1

            i += 1

    def export_struct_without_pseudoknots(self):
        """Returns the structure string with all pseudoknot base pairs removed.
        Use this for export to structural feature discovery.

        :rtype: string
        :returns: The secondary structure in dot-paren notation with all
            pseudoknotted base pairs removed.
        """
        struct = ''.split(self.secstruc)
        for i, char in enumerate(struct):
            if char in Constants.non_pknot_chars:
                struct[i] = '.'
        return ''.join(struct)

    def export_struct_pseudoknots(self):
        """Opposite of ``export_struct_without_pseudoknots()``. Export only
        pseudoknotted base pairs.

        :rtype: string
        :returns: The secondary structure in dot-paren notation where only
            pseudoknotted base pairs are retained.
        """
        struct = ''.split(self.secstruc)
        for i, char in enumerate(struct):
            if char in ['(', ')']:
                struct[i] = '.'
        return ''.join(struct)

    def base_pairs_by_pseudoknot_layer(self):
        """Returns a dict of lists of base pairs, one for each level of
        pseudoknottedness.

        Base pairs are indexed from 0.

        >>> struct = PseudoknotSecstruc('(([[))]]')
        >>> struct.base_pairs_by_pseudoknot_layer()
        {')': [(0, 5), (1, 4)], ']': [(2, 7), (3, 6)]}

        :rtype: dict(list(tuple(int, int)))
        :returns: A dictionary of base pair lists, where in each list is
            a set of base pairs at the same level of pseudoknottedness.
            The dictionary keys are **closing symbols** from the dot-paren
            representation used in this structure.

        :raises: ValueError
        """
        bp_at_levels = {')': []}  # For each closing character, a list
                                     # of closed base pairs.

        bp_stacks = {')': []}  # For each closing character, a stack of
                                  # currently open base pairs.
        levels = ['(']  # The opening characters for each encountered pknot
                        # level.

        for i, char in enumerate(self.secstruc):
            if char == '.':
                continue
            elif char in Constants.open_symbols:
                if char in levels:
                    bp_stacks[Constants.symbols[char]].append(i)
                else:
                    end_char = Constants.symbols[char]
                    bp_stacks[end_char] = [i]
                    bp_at_levels[end_char] = []
                    levels.append(char)
            elif char in Constants.close_symbols:
                if char in bp_stacks:
                    bp_at_levels[char].append((bp_stacks[char].pop(), i))
                else:
                    raise ValueError('Found closing symbol %s for unopened pair at position %d!' % (char,i))
            else:
                raise ValueError('Encountered invalid structure character %s at %i!' % (char, i))

        open_levels = [ l for l in levels if bp_stacks[Constants.symbols[l]] != [] ]
        if open_levels != []:
            raise ValueError('Invalid structure: not all base pairs have been closed (have you filetered tangled regions?), open levels: %s' % str(open_levels))

        # Reverse the stacks
        bp_at_levels_reversed = {
            level: list(reversed(bp_at_levels[level])) for level in bp_at_levels
        }

        return bp_at_levels_reversed

    def get_pseudoknot_indices(self):
        """Returns a list of base pairs that participate in pseudoknots in the
        structure.

        It's just a wrapper for :func:`get_pseudoknot_pairs`.

        .. warning::

          May be deprecated in the (near) future.

        :type bplist: list(tuple(int, int))
        :param bplist: A list of base pairs.

        :rtype: list(tuple(int, int))
        :returns: A list of base pairs such that they participate in a pseudoknot.
        """
        bplist = [bp for bp in self.get_base_pair_indices()]
        pknot = get_pseudoknot_pairs(bplist)
        return pknot


###############################################################################

#
# Some functions for manipulating secondary structures.
#


def contains_pseudoknot(bplist):
    """Checks if a list of base pairs contains a pseudoknot.

    :type bplist: list(tuple(int, int))
    :param bplist: A list of base pairs.

    :rtype: Boolean
    :returns: True if the base pairs in the list are not well-nested (the
        secondary structure represented by the base pairs contains a pseudoknot).
    """
    for a1, a2 in bplist:
        for b1, b2 in bplist:
            if (a1 < b1 < a2 < b2) \
                    or (a2 < b1 < a1 < b2) \
                    or (a2 < b2 < a1 < b1):
                return True
    return False
    
    
def get_pseudoknot_pairs(bplist):
    """Returns a list of base pairs that participate in pseudoknots.

    :type bplist: list(tuple(int, int))
    :param bplist: A list of base pairs.

    :rtype: list(tuple(int, int))
    :returns: A list of base pairs such that they participate in a pseudoknot.
    """
    result = []
    for a1, a2 in bplist:
        for b1, b2 in bplist:
            if (a1 < b1 < a2 < b2) or (a2 < b1 < a1 < b2) \
                    or (a2 < b2 < a1 < b1) or (b1 < a1 < b2 < a2) :
                if (a1, a2) not in result:
                    result.append((a1, a2))
    return result
                

def make_dotbracket_from_bplist(length, bplist, brackets=('(', ')'),
                                no_pknots=False):
    """Takes a structure as a list of basepairs and converts it to
    dot-paren notation.

     >>> bplist = [(2,7),(4,5),(10,15),(17,20)]
     >>> make_dotbracket_from_bplist(bplist)
     ((...).(((.((....))).)))

     Pseudoknotted structures raise a :class:`PseudoknotError` if ``no_pknots``
     is set.

    :type length: int
    :param length: How long the dot-bracket output should be (there may be
        unpaired positions *after* the last base pair closes).

    :type bplist: list(tuple(int, int))
    :param bplist: The secondary structure represented as a list of base pairs.
        The bplist indices are indices **into the resulting string**, indexed
        from **0**. If the structure contains a pseudoknot, a
        ``PseudoknotError`` is raised.

    :type brackets: tuple(char, char)
    :param brackets: If ``no_pknots`` is set, will use this tuple of left and
        right bracket to mark base pairs in the dot-paren string. If pseudoknots
        are permitted, has no effect.

    :type no_pknots: bool
    :param no_pknots: If set, will raise a :class:`PseudoknotError` when the
        base pair list contains a pseudoknot.

    :rtype: string
    :returns: The given structure in dot-paren notation.

    :raises: PseudoknotError
    """
    if no_pknots:
        if contains_pseudoknot(bplist):
            raise PseudoknotError("Base pairs contain a pseudoknot: "
                                  "%s" % str(bplist))
        secstruc = ['.']*length
        for bp in bplist:
            secstruc[bp[0]] = brackets[0]
            secstruc[bp[1]] = brackets[1]
        return ''.join(secstruc)

    else:

        pkr = PKnotResolver()
        bp_coloring = pkr.color_bplist(bplist)
        idx2bracket = { i : color
                        for color in bp_coloring
                        for bp in bp_coloring[color]
                        for i in bp }
        idx2bracket = {}
        for color in bp_coloring:
            complementary_bracket = Constants.symbols[color]
            for bp in bp_coloring[color]:
                idx2bracket[bp[0]] = color
                idx2bracket[bp[1]] = complementary_bracket

        #logging.debug('Received coloring: %s' % str(bp_coloring))
        #logging.debug('idx2bracket: %s' % str(idx2bracket))

        secstr = ['.' for _ in xrange(length)]
        for i in xrange(length):
            if i in idx2bracket:
                secstr[i] = idx2bracket[i]

        secstr_string = ''.join(secstr)

        #logging.debug('Output secstr: %s' % secstr_string)
        return secstr_string


def make_bplist_from_dotbracket(secstruc):
    """Reads a secondary structure object and returns a list of base pair
    indices. Can deal with pseudoknots.

    The following are equivalent:

    >>> make_bplist_from_dotbracket(secstruc)
    >>> secstruc.get_base_pair_indices_list()

    :type secstruc: Secstruc
    :param secstruc: The secondary structure object from which we wish
        to obtain a base pair list. Can also deal with a dot-paren string.

    :rtype: list(tuple(int, int))
    :returns: A list of the base pairs in the given structure.
    """
    result = []
    if isinstance(secstruc, Secstruc):
        return secstruc.base_pairs

    for bp in PseudoknotSecstruc(secstruc).get_base_pair_indices():
        result.append((bp[0]+1, bp[1]+1))
    return result


def has_overlapping_basepairs(bplist):
    """Checks whether there are base pairs in the list that share a residue.
    This should not happen, at least not in structures processable by this
    module.

    :type bplist: list(tuple(int, int))
    :param bplist: A secondary structure represented as a list of base pairs.

    :rtype: Boolean or int
    :returns: If the structure has no overlapping base pairs, returns False.
        If such a base pair is found, returns the index of its overlapping
        member (if there are duplicates, returns the start of the bae pair).
    """
    indices = {}
    for bp in bplist:
        if indices.has_key(bp[0]): return bp[0]
        if indices.has_key(bp[1]): return bp[1]
        indices[bp[0]]=True
        indices[bp[1]]=True
    return False
    

def filter_by_mask(structure, mask, keep_outpairs=False):
    """Filters all base pairs that have ends outside the mask.

    Filters all base pairs that have ends outisde the mask. Outpairs (pairs
    with one masked and one un-masked end) are kept if ``keep_outpairs`` is set
    and discarded otherwise.

    Returns a secondary structure with only the correct BPs left.

    :type structure: Secstruc
    :param structure: The original structure from which to filter.

    :type mask: list(bool)
    :param mask: A mask (Boolean) array of the same length as the structure.
        Only pairs with both (or at least one, if ``keep_outpairs`` is set) ends
        on a ``True`` index of the mask are retained.

    :type keep_outpairs: bool
    :param keep_outpairs: If set, will not discard base pairs with one end
        inside the mask and one outside.

    :rtype: Secstruc
    :returns: A secondary structure with base pairs outside the mask filtered
        out. The type of the structure is retained (pseudoknotted structures
        stay pseudoknotted, non-pknot structures stay non-pknot).
    """
    if not isinstance(structure, Secstruc):
        raise TypeError('Supplied structure not a Secstruc object:'
                        '%s' % str(type(structure)))

    bp_list = structure.base_pairs
    retain_bplist = []
    for bp in bp_list:
        #logging.debug('  At base pair: %s' % str(bp))
        if mask[bp[0]] and mask[bp[1]]:
            retain_bplist.append(bp)
        elif mask[bp[0]] != mask[bp[1]] and keep_outpairs:
            retain_bplist.append(bp)

    logging.debug('Creating dot-bracket from bplist %s' % str(retain_bplist))
    dotbracket = make_dotbracket_from_bplist(len(structure), retain_bplist)
    output_structure = structure.__class__(dotbracket, structure.indices)

    #logging.debug('Output structure type: %s' % type(output_structure))
    return output_structure


def validate_secstruc_string(secstruc, raise_exception=False):
    structure = PseudoknotSecstruc(secstruc)

    dpchars = set(secstruc)
    bracket_levels = [(b, Constants.symbols[b])
                      for b in Constants.bracket_priorities if b in dpchars]

    is_valid = True
    invalid_levels = []
    for brackets in bracket_levels:
        if not structure.validate(diagnostics=False, brackets=brackets):
            invalid_levels.append(brackets)
            is_valid = False

    if not is_valid:
        if raise_exception:
            raise ValueError('Validating secstruc string: INVALID\t'
                             'Levels with invalid BPs: {0}'.format(invalid_levels))
        else:
            logging.error('Validating secstruc string: INVALID\t'
                          'Levels with invalid BPs: {0}'.format(invalid_levels))

    return is_valid
