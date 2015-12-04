"""
This module contains utility functions that manipulate alignment data and
sequence/structure data with respect to an alignment: expanding/contracting
according to where gaps in the alignment are.
"""
import logging
import os
import random
import string
from Bio import AlignIO
from Bio import SeqIO
import Bio
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Phylo.TreeConstruction import DistanceCalculator
import operator
import itertools
import collections
from Bio.Seq import Seq
import copy
from Bio.SeqRecord import SeqRecord
import Secstruc as SS
#from rPredictor import ALLOWED_BPS, GAP_CHARS
# import matplotlib.pyplot as plt


__author__ = "Jan Hajic jr."


def seq_aln2seq(aln_seq_string):
    """Given a sequence string, removes all gaps.

    :type aln_seq_string: str
    :param aln_seq_string: The aligned string from which to remove gaps.
        Gap characters are ``-`` and ``_`` (rare, but happens in some tools).

    :rtype: string
    :returns: The resulting sequence string.
    """
    ALN_GAP = ['-', '_']

    # if isinstance(aln_seq_rec, SeqRecord)... etc.
    # aln_seq_string = str(aln_seq_rec.seq)

    seq_chars = []
    for c in aln_seq_string:
        if c not in ALN_GAP:
            seq_chars.append(c)

    seq_string = ''.join(seq_chars)
    return seq_string


def struct_seq2aln(struct, aln_seq):
    """Given a structure, and an aligned sequence,
    expands the structure dot-paren string to match the alignment.
    Gaps in alignment will be represented as dots in the aligned
    structure.

    :type struct: str
    :param struct: A dot-paren representation of a secondary structure. It
        should correspond to the ``aln_seq`` sequence.

    :type aln_seq: str
    :param aln_seq: An aligned sequence. Should be of the same length as
        the ``struct`` dot-paren string.

    :rtype: str
    :returns: The resulting aligned structure dot-paren string.
    """

    # Why was this even here? More importantly, why did it work?
    #assert(len(struct) == len(aln_seq))

    ALN_SEQ_TO_STRUCT_MAP = { '-':'.', '_':'.' }
    ALN_GAP = ALN_SEQ_TO_STRUCT_MAP.keys()

    # Extract dot-paren string from structure.
    structure_string = str(struct)

    # Expand dot-paren structure string
    output_struct_chars = []

    i = 0 # I is the index into the original sequence (no gaps), used to
          # access structure
    for j in xrange(len(aln_seq)): # J is the index into the aligned sequence
        if aln_seq[j] in ALN_GAP: # Gap in alignment - need to expand structure
            aln_struct_gap_char = ALN_SEQ_TO_STRUCT_MAP[aln_seq[j]]
            output_struct_chars.append(aln_struct_gap_char)
        else:
            output_struct_chars.append(structure_string[i])
            i += 1

    output_struct_string = ''.join(output_struct_chars)
    return output_struct_string


def struct_aln2seq(aln_struct, aln_seq, no_check_canonical=False):
    """Given a structure and an aligned sequence, contracts the structure
    so that base pairs corresponding to gaps in the given sequence
    alignment are removed and positions in the structure dot-paren
    string corresponding to gaps are cut out.

    This function checks the elementary problems of fitting structures
    that come from different molecules via an alignment to the target
    sequences:

    * base pairs of the original structure that have one or both ends in a gap
      of the target structure (while in the aligned *source* structure, they
      pointed to nucleotides). The following aligned source structure and
      sequence::

        AACG--U-GGC-GACUU
        ((((..(..)..).)))

      with the target sequence aligned to the source structure::

        CACG-GU-G-CAGAC-U
        ((((..(..)..).)))

      has this obvious problem at the position in the ``-G-`` and ``C-U``
      regions: a source base pair corresponding to a gap in the aligned target
      sequence. This function will filter out the participating base pairs::

        CACG-GU-G-CAGAC-U
        (.((........).).)

    * non-canonical base pairs can also form (see the outermost pair in the
      previous example). Unless the ``no_check_canonical`` flag is given,
      non-canonical base pairs will also be removed::

        CACG-GU-G-CAGAC-U
        ..((........).)..

    Finally, the contracted sequence/structure pair will be::

      CACGGUGCAGACU
      ..((.....).).

    and the return value::

      ..((.....).).

    :type aln_struct: str
    :param aln_struct: A dot-paren representation of a secondary structure
        of the sequence ``aln_seq``. The structure should be contracted so that
        positions corresponding to gaps in the alignment from which ``aln_seq``
        is taken are removed.

        .. note::

          The structure at this point might NOT be a "real" structure for the
          given sequence - this function makes sure that the *output* is a
          valid dot-paren representation of the ``aln_seq``'s secondary
          structure, but the input structure may come from a different molecule.

    :type aln_seq: str
    :param aln_seq: The aligned sequence according to which the ``aln_struct``
        dot-paren string will be contracted.

    :type no_check_canonical: Boolean
    :param no_check_canonical: If this flag is set to True, will not remove
        all non-canonical base pairs from the structure. (Default: ``False``)

    :rtype: str
    :returns: The dis-aligned structure dot-paren string.
    """

    # Sanity checks
    if len(aln_struct) != len(aln_seq):
        raise ValueError("Different length of structure (%d) and sequence (%d)!" % (len(aln_struct), len(aln_seq)))

    SEQ_NOPAIR_CHAR = '.'
    ALN_GAP_SEQ_TO_STRUCT_MAP = { '-': '.', '_': '.' }

    aln_struct_chars = list(str(aln_struct))

    # DEBUG
    #aln_seq_chars = list(aln_seq)
    #for i in xrange(len(aln_struct_chars)):
    #    print aln_seq_chars[i], aln_struct_chars[i]
    #print no_check_canonical, "NO_CHECK_CANONICAL on, will not filter out non-canonical base pairs."

    # Process base pairs
    logging.debug('Computing aligned structure bplist...')
    bplist = SS.make_bplist_from_dotbracket(aln_struct)

    logging.debug('Iterating through base pairs.')
    for first, second in bplist:
        # Correct for base pair indices starting from 1
        first -= 1
        second -= 1

        if aln_seq[first] in ALN_GAP_SEQ_TO_STRUCT_MAP:
            # DEBUG
            #print "Removing base pair %d-%d, where %d is a gap (%s / %s) and %d is %s / %s" % (first, second, first, aln_seq[first], aln_struct_chars[first], second, aln_seq[second], aln_struct_chars[second])
            aln_struct_chars[first] = ALN_GAP_SEQ_TO_STRUCT_MAP[aln_seq[first]]
            aln_struct_chars[second] = ALN_GAP_SEQ_TO_STRUCT_MAP[aln_seq[first]]
        elif aln_seq[second] in ALN_GAP_SEQ_TO_STRUCT_MAP:
            # DEBUG
            #print "Removing base pair %d-%d, where %d is a gap (%s / %s) and %d is %s / %s" % (first, second, second, aln_seq[second], aln_struct_chars[second], first, aln_seq[first], aln_struct_chars[first])
            aln_struct_chars[first] = ALN_GAP_SEQ_TO_STRUCT_MAP[aln_seq[second]]
            aln_struct_chars[second] = ALN_GAP_SEQ_TO_STRUCT_MAP[aln_seq[second]]
        elif (no_check_canonical is False) and ((aln_seq[first], aln_seq[second]) not in ALLOWED_BPS):
            # DEBUG
            #print "Removing base pair %d-%d, which is not allowed (%d is %s, %d is %s)." % (first, second, first, aln_seq[first], second, aln_seq[second])
            aln_struct_chars[first] = SEQ_NOPAIR_CHAR
            aln_struct_chars[second] = SEQ_NOPAIR_CHAR

    # Cut out gap positions
    struct_chars = []
    for i, nt in enumerate(aln_seq):
        if nt not in ALN_GAP_SEQ_TO_STRUCT_MAP:
            struct_chars.append(aln_struct_chars[i])

    struct_string = ''.join(struct_chars)
    return struct_string


def contract_regions_wrt_aln(regions, aln_seq):
    """Given a list of regions and an alignment sequence from which the regions
    were taken, returns the list of regions recomputed for indices that are
    valid after contracting the sequence back to its unaligned form.

    :type regions: list(tuple(int, int))
    :param regions: A list of regions to contract.

    :type aln_seq: str
    :param aln_seq: The aligned sequence into which the input regions are valid.

    :rtype: list(tuple(int, int))
    :returns: A list of the input regions with their ``(start, end)`` indices
        recomputed to be valid for ``aln_seq`` after contracting it to the
        unaligned version.
    """
    logging.debug('Contracting regions wrt alignment.')
    logging.debug('Aligned sequence: {0}'.format(aln_seq))
    logging.debug('Regions to contract: {0}'.format(regions))

    GAP_CHARS = set(['-', '_'])

    # ...so that we do not unnecessarily generate gap counts.
    if len(regions) == 0:
        return []

    # Prepare gap counts
    logging.debug('Preparing gap counts...')
    gap_counter = 0
    gap_counts = []
    for char in aln_seq:
        if char in GAP_CHARS:
            gap_counter += 1
        gap_counts.append(gap_counter)

    # Process regions
    logging.debug('Processing regions (total: %d)...' % len(regions))
    output_regions = []
    for region in regions:
        start, stop = region

        if aln_seq[start] in GAP_CHARS:
            while aln_seq[start] in GAP_CHARS and start < len(gap_counts) - 1:
                start += 1
            # Problem: what if there are no non-gaps up to the end? We get
            # start = len(gap_counts) - 1, which is aln_seq[-1], and stop.
            # If this is also gap, the region doesn't map to anything in the
            # uncontracted sequence.
            if aln_seq[start] in GAP_CHARS:
                continue
        c_start = start - gap_counts[start]

        if aln_seq[stop] in GAP_CHARS:
            while aln_seq[stop] in GAP_CHARS and stop > 0:
                 stop -= 1
            # Same problem here: what if region falls into a gaps-only start
            # region of the aln seq?
            if aln_seq[stop] in GAP_CHARS:
                continue
        c_stop = stop - gap_counts[stop]

        logging.debug('From (%d, %d) to (%d, %d) over %s' % (start, stop,
                                                             c_start, c_stop,
                                                             aln_seq[start:stop+1]))
        if c_start <= c_stop:
            output_regions.append((c_start, c_stop))

    return output_regions


def expand_regions_wrt_aln(regions, aln_seq):
    """Complementary to :func:`contract_regions_wrt_aln`. Given a set of
    regions that point to the unaligned equivalent of ``aln_seq``, expands the
    region indices so that they correspond to indices to the aligned sequence.

    :type regions: list(tuple(int, int))
    :param regions: The set of regions to be expanded.

    :type aln_seq: string
    :param aln_seq: The aligned form of the sequence for which the ``regions``
        are valid.

    :rtype: list(tuple(int, int))
    :returns: A list of the expanded regions.
    """
    GAP_CHARS = set(['-', '_'])

    output_regions = []
    if len(regions) == 0:
        return output_regions

    gap_count = 0

    regions_iter = regions.__iter__()
    unaligned_idx = 0
    current_region = regions_iter.next()
    output_region = [current_region[0], current_region[1]]

    for aligned_idx, a in enumerate(aln_seq):
        if a in GAP_CHARS:
            gap_count += 1
            continue
        if unaligned_idx == current_region[0]:
            output_region[0] = aligned_idx
        if unaligned_idx == current_region[1]:
            output_region[1] = aligned_idx
            output_regions.append(tuple(output_region))
            try:
                current_region = regions_iter.next()
                output_region = [current_region[0], current_region[1]]
            except StopIteration:
                break
        unaligned_idx += 1

    return output_regions


def recode_regions(regions, source_aln_seq, tgt_aln_seq):
    """Given a set of regions for a source sequence, its aligned form
    and a target sequence aligned to the same form, re-codes the given
    regions into the target unaligned sequence.

    >>> regions = [(1,3), (6,6), (9,13)]
    >>> source_aln_seq = 'AAC-GG-UCAAGUG--UCGA-U---'
    >>> tgt_aln_seq =    'AA-AGCCU--AG-GACUCG-AUUGC'
    >>> recode_regions(regions, source_aln_seq, tgt_aln_seq)
    [(1,3), (8, 13)]

    Why is this the result? Let's go step by step. In the unaligned form of
    ``src_aln_seq``, we have these regions (denoted by dashes)::

      AAACGGUCAAGUGUCGAU
      .---..-..------...

    In the aligned form, they turn to::

      AAC-GG-UCAAGUG--UCGA-U---
      .----...-..-------.......

    Applied to the target aligned sequence, they are::

      AA-AGCCU--AG-GACUCG-AUUGC
      .----...-..-------.......

    And after contracting, we obtain::

      AAAGCCUAGGACUCGAUUGC
      .---....------......

    and thus ``[(1,3), (8, 13)]``. Notice the disappearing region that
    was mapped over to gaps only.


    :type regions: list(tuple(int,int))
    :param regions: A list of regions valid for the unaligned form of
        ``source_aln_seq``.

    :type source_aln_seq: string
    :param source_aln_seq: The aligned form of the source sequence. The
        given ``regions`` are valid for its unaligned form.

    :type tgt_aln_seq: string
    :param tgt_aln_seq: The aligned form of the target sequence. The ``regions``
        will be re-coded to be valid for the unaligned form of this sequence.

    :rtype: list(tuple(int, int))
    :returns: A list of regions valid for the unaligned form of ``tgt_aln_seq``.
    """
    assert(len(source_aln_seq) == len(tgt_aln_seq))

    expanded_regions = expand_regions_wrt_aln(regions, source_aln_seq)
    recoded_regions = contract_regions_wrt_aln(expanded_regions, tgt_aln_seq)

    return recoded_regions


def expand_structure_wrt_aln(structure, aln_seq):
    """Given a Secstruc object, returns another with its indices and
    secstr expanded according to the aligned (i.e. with gaps) sequence.

    Gaps will correspond to dots in the dot-paren representation.

    The situation with the ``indices`` is more complicated. The gap positions
    will be set **to the value of -1 in the indices array**,
    to reflect that the gap positions are **not** a part of the molecule.
    Indices provide the residue numbering; a residue's number does not change
    when placeholders are inserted around it.

    >>> structure = SS.Secstruc('(((...))).', range(10))
    >>> aln_seq = 'AAC--GUUG-UU--C'
    >>> newstruct = expand_structure_wrt_aln(structure, aln_seq)
    >>> newstruct.secstruc
    '(((.....).))...'
    >>> newstruct.indices
    [0, 1, 2, -1, -1, 3, 4, 5, 6, -1, 7, 8, -1, -1, 9]

    :type structure: rPredictor.Secstruc.Secstruc
    :param structure: The secondary structure to expand.

    :type aln_seq: str
    :param aln_seq: The (gapped) sequence according to which the structure
        should be expanded.

    :rtype: rPredictor.Secstruc.Secstruc
    :returns: The expanded secstruc.
    """
    dpstring = structure.secstruc
    indices = structure.indices

    expanded_dpstring = struct_seq2aln(dpstring, aln_seq)
    expanded_indices = [ -1 for _ in aln_seq ]
    original_position = 0
    for expanded_position, char in enumerate(aln_seq):
        if char not in GAP_CHARS:
            expanded_indices[expanded_position] = indices[original_position]
            original_position += 1

    # Secstruc or PseudoknotSecstruc
    expanded_structure = structure.__class__(expanded_dpstring,
                                             expanded_indices)
    return expanded_structure


def contract_alignment(alignment):
    """Creates a new alignment that discards gap-only positions of the given
    alignment. Intended for situations where a sub-alignment has been
    constructed and some positions have been left without nucleotides.
    """
    if len(alignment) == 0:
        return copy.deepcopy(alignment)

    new_seqs_chars = [[] for _ in alignment]

    total_discarded = 0
    total_retained = 0
    for position in xrange(len(alignment[0])):
        column = alignment[:,position]
        has_non_gap = False
        for c in column:
            if c not in GAP_CHARS:
                has_non_gap = True
        if has_non_gap is True:
            total_retained += 1
            for i, c in enumerate(column):
                new_seqs_chars[i].append(c)
        else:
            total_discarded += 1

    new_seq_strings = [''.join(seq_chars) for seq_chars in new_seqs_chars]
    new_seqs = [Seq(s, SingleLetterAlphabet())
                for s in new_seq_strings]

    new_records = []
    for record, seq in itertools.izip(alignment, new_seqs):
        new_record = SeqRecord(seq,
                               id=record.id,
                               name=record.name,
                               description=record.description,
                               dbxrefs=record.dbxrefs,
                               features=record.features,
                               annotations=record.annotations,
                               letter_annotations=record.letter_annotations)
        new_records.append(new_record)

    output_aln = MultipleSeqAlignment(new_records)
    return output_aln


###############################################################################

# Retrieving things from alignments


def find_seq_record_in_aln(aln_seq_records, name):
    """Returns the sequence record for the given name
    from the given alignment.

    If sequence record with given name is not found, returns ``None``.

    :type aln_seq_records: list(Bio.Seq.SeqRecord)
    :param aln_seq_records: The list of aligned sequence records obtained
        from an alignment.

    :type name: string
    :param name: The ``id`` member of the sequence we want to find. (This
        is, for example, the name under which it lives in a CLUSTAL alignment
        file.)

    :rtype: Bio.Seq.SeqRecord
    :returns: The sequence record with the given name, or ``None`` if such
        record is not found in the ``aln_seq_records``.
    """
    for aln_seq_record in aln_seq_records:
        if aln_seq_record.id == name:
            return aln_seq_record

    return None


def extract_seq_string_from_aln(aln_seq_records, name):
    """Extracts sequence string for the given name from
    the given alignment records (or alignment).

    If sequence record with given name is not found, raises ``ValueError``.

    :type aln_seq_records: list(Bio.Seq.SeqRecord)
    :param aln_seq_records: The list of aligned sequence records obtained
        from an alignment.

    :type name: string
    :param name: The ``id`` member of the sequence we want to find. (This
        is, for example, the name under which it lives in a CLUSTAL alignment
        file.)

    :rtype: string
    :returns: The sequence with the given name.

    :raises: ValueError
    """
    if isinstance(aln_seq_records, Bio.Align.MultipleSeqAlignment):
        aln_seq_records = [ a for a in aln_seq_records ]

    aln_seq_record = find_seq_record_in_aln(aln_seq_records, name)
    if not aln_seq_record:
        raise ValueError("Sequence with name %s not found in alignment. (Aln. names: %s)" % (name, str([rec.id for rec in aln_seq_records])))

    return str(aln_seq_record.seq)


def read_alignment_records(alignment_file, fmt="clustal"):
    """Returns a list of seq records contained in the alignment file.
    Assumes the file only contains one alignment. (Stick with that anyway.)
    Assumes CLUSTAL alignment format by default.

    :type alignment_file: string
    :param alignment_file: The file where the alignment is stored.

    :type fmt: string
    :param fmt: A format recognized by ``Bio.AlignIO`` corresponding to the
        format of the ``alignment_file``.

    :rtype: list(Bio.Seq.SeqRecord)
    :returns: A list of Biopython sequence records in the alignment.
    """
    with open(alignment_file) as aln_handle:
        aln_generator = AlignIO.parse(aln_handle, fmt)
        aln = aln_generator.next()

    aln_seq_records = []
    for aln_seq_record in aln:
        aln_seq_records.append(aln_seq_record)


    return aln_seq_records


def read_alignment(alignment_file, fmt="clustal"):
    """Returns alignment from file. Assumes the file only contains one
    alignment (stick with that, anyway). Assumes CLUSTAL format by default.

    :type alignment_file: string
    :param alignment_file: The file where the alignment is stored.

    :type fmt: string
    :param fmt: A format recognized by ``Bio.AlignIO`` corresponding to the
        format of the ``alignment_file``.

    :rtype: Bio.Align.MultipleSeqAlignment

    """
    with open(alignment_file, 'r') as aln_handle:
        aln_gen = AlignIO.parse(aln_handle, fmt)
        aln = aln_gen.next()

    aln._alphabet = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.IUPACAmbiguousRNA)

    return aln


def random_k_from_alignment(alignment, k=10, contract=False):
    """Return a new alignment with K records chosen randomly from the given
    alignment. If the contract flag is set, contracts the output alignment."""
    if k > len(alignment):
        logging.warn('Trying to sample {0} records from alignment with only '
                     '{1} records. Retrieving only {1} records.'
                     ''.format(k, len(alignment)))
        k = len(alignment)

    indices = range(len(alignment))
    random.shuffle(indices)

    use_indices = indices[:k]
    records = [alignment[i] for i in use_indices]
    output = MultipleSeqAlignment(records)

    if contract is True:
        output = contract_alignment(output)

    return output


def read_records_from_aln_or_seqs(filename, fmt):
    """If 'fmt' is clustal, reads records from alignment file 'filename'.
    If 'fmt' is fasta, reads the records through SeqIO."""
    raise NotImplementedError()

###############################################################################

# Dealing with duplicate IDs. Helper functions for aligning.


def find_duplicate_records(*records):
    """Returns a dictionary grouped by ID of records that have the same
    ID."""
    records_by_ids = collections.defaultdict(list)
    ids_with_duplicates = set()
    for rec in itertools.chain(*records):
        if rec.id in records_by_ids:
            ids_with_duplicates.add(rec.id)
        records_by_ids[rec.id].append(rec)
    records_with_duplicates = {recid: records_by_ids[recid]
                               for recid in ids_with_duplicates}
    return records_with_duplicates


def deduplicate_records(*records):
    """Renames sequence records with duplicate IDs and returs a list of records
    with those re-named IDs."""
    output_records = []
    id_counts = {}
    for rec in itertools.chain(*records):
        if rec.id not in id_counts:
            output_records.append(rec)
            id_counts[rec.id] = 1
        else:
            new_id = rec.id + '.dupl-{0}'.format(id_counts[rec.id] - 1)
            new_rec = copy.deepcopy(rec)
            new_rec.id = new_id
            output_records.append(new_rec)
            id_counts[rec.id] += 1
    return output_records


def are_duplicates_in_seqfiles(files, formats):
    """Checks whether there are records with duplicate IDs across the given
    files. ``formats`` is a list that says which format the file is in --
    accepts fasta for sequences, clustal for alignments."""
    logging.info('Checking for duplicate IDs (clustalw2 safeguard)...')
    all_records = []
    for f, fmt in zip(files, formats):
        if fmt == 'fasta':
            records = list(SeqIO.parse(open(f), fmt))
        elif fmt == 'clustal':
            records = [r for r in read_alignment(f, fmt)]
        else:
            raise ValueError('Unsupported format: {0}!'.format(fmt))
        all_records.append(records)
    duplicates = find_duplicate_records(*all_records)
    if len(duplicates) == 0:
        logging.info('No duplicate records found.')
        return False
    else:
        log_string = 'ClustalW2 error prevention found the following ' \
                     'duplicate records: '
        for d_id in duplicates:
            log_string += '\n  {0}: {1} x'.format(d_id, len(duplicates[d_id]))
        logging.warn(log_string)
        return True


def rename_duplicates_from_file(input_file, other_files_to_check,
                                input_fmt, other_fmts,
                                output_file=None, output_fmt=None):
    """Rewrites the input sequence records file to the output so that the
    duplicate records in the input file are renamed to not cause conflict with
    records in other_files_to_check. Returns the set of records from the input
    file so that all records that have duplicates in ``other_files_to_check``
    are renamed to avoid collisions with the other files.

    If ``output_file`` and ``output_fmt`` are specified, outputs the renamed
    records to the given file in the given format.

    The formats are either 'clustal' or 'fasta'. There must be exactly as many
    'other_fmts' as there are 'other_files_to_check'.
    """
    if len(other_files_to_check) != len(other_fmts):
        raise ValueError('{0} files to check but {1} formats!'
                         ''.format(len(other_files_to_check), len(other_fmts)))
    input_records = read_records_from_aln_or_seqs(input_file, input_fmt)
    check_records = [read_records_from_aln_or_seqs(f, fmt)
                     for f, fmt in zip(other_files_to_check, other_fmts)]

    duplicates = find_duplicate_records(input_records, check_records)
    duplicates_check = find_duplicate_records(check_records)
    duplicates_input = find_duplicate_records(input_records)

    # Duplicates in check-records: only rename those in the input records.
    # - renaming strategy??

    # Output?

    # Return renamed.

###############################################################################

# Aligning


def align_to_profile(profile_alignment_file, sequence_fst, tmp_dir=''):
    """Aligns the sequence to the given profile alignment file. Returns
    the complete alignment.

    :type profile_alignment_file: str
    :param profile_alignment_file: The file containing the profile alignment
        to which the input sequence should be aligned. CLUSTAL format is
        assumed.

    :type sequence_fst: str
    :param sequence_fst: The FASTA file from which to read the sequence for
        which we want to select a template.

        Make sure the FASTA header in this file cannot match the template
        names in ``template_aln_file``.

    :type tmp_dir: str
    :param tmp_dir: A directory for temporary files where the temporary
        alignment created by aligning the input sequence to the template
        selection alignment is created. (The temporary file does not
        "survive" this method.)

    :rtype: Bio.Align.MultipleSeqAlignment
    :returns: The alignment object with the input sequence aligned to the
        profile.
    """
    # Run clustalw2 profile alignment,
    # parse output as Alignment object

    random_string = ''.join([random.choice(string.ascii_lowercase + string.digits) for _ in xrange(16)])
    tmp_aln_name = os.path.join(tmp_dir, 'cp_predict.' + random_string + '.tmp.aln')

    logging.debug('Creating tmp alignment file: {0}, using profile {1}'
                  ''.format(tmp_aln_name, profile_alignment_file))

    are_duplicates_in_seqfiles([profile_alignment_file, sequence_fst],
                               ['clustal', 'fasta'])

    # Find duplicates from alignment file, create temp aln file. Seqfile
    # sequences will *not* be renamed. (NOT IMPLEMENTED)

    # Check if clustal exists
    CLUSTAL_CMD = rPredictor.get_clustal_command()
    logging.info('Will be running clustal: {0}'.format(CLUSTAL_CMD))
    cline = ClustalwCommandline(CLUSTAL_CMD,
                                profile2=sequence_fst,
                                profile1=profile_alignment_file,
                                sequences=True,
                                outfile=tmp_aln_name)
    cline()
    alignment = read_alignment(tmp_aln_name)

    # Rename duplicates back to original (this may fail, have to check).
    # Note that this may produce unusable TTAs... maybe don't restore the names?

    os.remove(tmp_aln_name)

    # clustalw2 also creates a guide tree, which we do not need.
    tmp_dnd_name = sequence_fst[:-4] + '.dnd'
    if os.path.isfile(tmp_dnd_name):
        os.remove(tmp_dnd_name)
    else:
        logging.warn('ClustalW2 did not create tmp guide tree file, something '
                     'strange is going on. Trying to continue.')

    return alignment


def align_msa(sequences_fst, tmp_dir=''):
    """Aligns the sequences from the given file. Returns
    the complete alignment.

    :type sequences_fst: str
    :param sequences_fst: The FASTA file from which to read the sequence which
        we want to align.

    :type tmp_dir: str
    :param tmp_dir: A directory for temporary files where the temporary
        alignment created by aligning the input sequence to the template
        selection alignment is created. (The temporary file does not
        "survive" this method.)

    :rtype: Bio.Align.MultipleSeqAlignment
    :returns: The alignment object with the input sequence aligned to the
        profile.
    """
    # Run clustalw2 multiple sequence alignment,
    # parse output as Alignment object

    random_string = ''.join([random.choice(string.ascii_lowercase + string.digits) for _ in xrange(16)])
    tmp_aln_name = os.path.join(tmp_dir, 'cp_predict.' + random_string + '.tmp.aln')

    logging.debug('Creating tmp alignment file: {0}, using sequences {1}'
                  ''.format(tmp_aln_name, sequences_fst))

    are_duplicates_in_seqfiles([sequences_fst],
                               ['fasta'])

    CLUSTAL_CMD = rPredictor.get_clustal_command()
    logging.info('Will be running clustal: {0}'.format(CLUSTAL_CMD))
    cline = ClustalwCommandline(CLUSTAL_CMD,
                                infile=sequences_fst,
                                align=True,
                                outfile=tmp_aln_name)
    cline()
    alignment = read_alignment(tmp_aln_name)

    os.remove(tmp_aln_name)

    # clustalw2 also creates a guide tree, which we do not need.
    tmp_dnd_name = sequences_fst[:-4] + '.dnd'
    if os.path.isfile(tmp_dnd_name):
        os.remove(tmp_dnd_name)
    else:
        logging.warn('ClustalW2 did not create tmp guide tree file {0}, '
                     'something strange is going on. Trying to continue.'
                     ''.format(tmp_dnd_name))

    return alignment


###############################################################################

# Distance-measuring tools


def distances_to_seq(alignment, sequence, distance_model='identity'):
    """A tool for computing not the complete sequence-sequence distance matrix,
    but only the distances to certain sequences.

    Beware: relies on a protected member of DistanceCalculator.

    :param alignment: A MultipleSeqAlignment object.

    :param sequence: A SeqRecord object. Must be of the same length as the
        records in the alignment.

    :param distance_model: One of either 'identity', 'blastn', or 'trans'.
        Defines the distance of a nucleotide pair. See
        Bio.Phylo.TreeConstruction.DistanceCalculator documentation.

    :returns: A list of distances between the given sequence and all sequences
        in the MSA, in the order in which the sequences are in the MSA.
    """
    dcalc = DistanceCalculator(distance_model)
    output = [dcalc._pairwise(sequence, msa_seq) for msa_seq in alignment]
    return output


def sort_by_distance(alignment, sequence, distance_model='identity'):
    """Sorts the sequence records in the alignment by distance from the given
    sequences. Returns the sorted list of (SeqRecord, distance) pairs.
    """
    distances = distances_to_seq(alignment, sequence,
                                 distance_model=distance_model)
    distance_record_pairs = [(record, distance)
                             for record, distance in zip(alignment, distances)]
    sorted_distance_record_pairs = sorted(distance_record_pairs,
                                          key=operator.itemgetter(1))
    return sorted_distance_record_pairs


def aln_from_closest_k(records, distances, k):
    """Given a set of sequence records and distances, creates an alignment
    from the first K closest sequence records. The distances are not assumed
    to be sorted, but the sequence records are assumed to all be of the same
    length, in order to create the alignment."""
    sorted_dr_pairs = sorted(zip(records, distances),
                             key=operator.itemgetter(1))
    sorted_records = map(operator.itemgetter(0), sorted_dr_pairs)
    closest_k_records = sorted_records[:k]
    #aln = MultipleSeqAlignment(closest_k_records)
    return closest_k_records


def aln_from_distance(records, distances, dist):
    """Given a set of sequence records and distances, creates an alignment from
    all records with distance less than the given number. The records in the
    alignment will NOT be sorted from closest to furthest.
    """
    return [r for r, d in itertools.izip(records, distances) if d < dist]


def aln_from_range(records, distances, min, max):
    """Given a set of sequence records and distances, creates an alignment from
    all records with distance d, such that ``min <= d < max``. The records in
    the alignment will NOT be sorted from closest to furthest.
    """
    return [r for r, d in itertools.izip(records, distances)
            if (min <= d) and (d < max)]


def distribute_by_distances(alignment, distances, bins, cumulative=False,
                            contract=True):
    """Retruns a list of multiple sequence alignments so that the distances
    of each of the i-th member of the list fall inside the range of (bins[i],
    bins[i+1]). If the ``cumulative`` flag is set, the i-th alignment will
    contain all the sequences from the previous bins.

    :param alignment: A MultipleSeqAlignment.

    :param distances: The distances (from some pivot sequence, not necessarily
        in the alignment). The i-th distance refers to the i-th record in the
        alignment.

    :param bins: A list of (k+1) bounds in ascending order. The j-th alignment
        returned will correspond to the interval between bins[j] and bins[j+1].

    :param cumulative: If set, the j-th alignment will also contain all
        sequences that fell into bins 0 ... j-1.

    :param contract: If set, will contract the output alignments so that there
        are no gap-only positions.

    :return: A list of MultipleSeqAlignment objects, one for each bin.
    """
    record_sets = []   # These alignments are *not* contracted, to be join-able
    records = [r for r in alignment]
    for min, max in zip(bins[:-1], bins[1:]):
        aln = aln_from_range(records, distances, min, max)
        record_sets.append(aln)
    # print 'Record sets: {0}'.format([len(r) for r in record_sets])
    if cumulative:
        cumulative_rsets = [record_sets[0]]
        # rint 'Cumulative rsets: {0}'.format(len(cumulative_rsets[0]))
        for i, rset in enumerate(record_sets[1:]):
            # print 'Rset {0}'.format(i)
            prev_rset = copy.deepcopy(cumulative_rsets[i])
            prev_rset.extend(rset)
            cumulative_rsets.append(prev_rset)
            # print 'Cumulative rsets: {0}'.format([len(r) for r in cumulative_rsets])
        record_sets = cumulative_rsets

    if contract is True:
        alns = [contract_alignment(MultipleSeqAlignment(rset))
                for rset in record_sets]
    else:
        alns = [MultipleSeqAlignment(rset) for rset in record_sets]

    logging.debug('Alignments: {0}'.format('\n'.join([str(a) for a in alns])))
    return alns


def distribute_by_closest(alignment, distances, sizes, contract=True):
    """Returns a list of Multiple Sequence Alignments so that each alignment
    consists of the K sequences with the lowest distance. For example, if the
    ``sizes`` parameter is ``[10, 20, 50, 100]``, the function will output
    a list of four MSAs, consisting of 10, 20, 50 and 100 closest sequences by
    distances. (The sizes do not have to be sorted.)

    >>> from Bio.Alphabet import generic_dna
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> a = SeqRecord(Seq("AAAA-CCCG-", generic_dna), id="Alpha")
    >>> b = SeqRecord(Seq("AAAA-CCGG-", generic_dna), id="Beta")
    >>> c = SeqRecord(Seq("A--ACCCGUU", generic_dna), id="Gamma")
    >>> align = MultipleSeqAlignment([a, b, c])
    >>> distances = [0.1, 0.15, 0.2]
    >>> alns = distribute_by_closest(align, distances, [1, 2, 3])
    >>> len(alns)
    3
    >>> [len(a) for a in alns]
    [1, 2, 3]
    >>> [[r.id for r in a] for a in alns]
    [['Alpha'], ['Alpha', 'Beta'], ['Alpha', 'Beta', 'Gamma']]
    >>> [len(a[0]) for a in alns]
    [8, 8, 10]

    :param alignment: A MultipleSeqAlignment.

    :param distances: The distances (from some pivot sequence, not necessarily
        in the alignment). The i-th distance refers to the i-th record in the
        alignment.

    :param sizes: A list of the output alignment sizes.

    :param contract: If set, will contract the output alignments so that there
        are no gap-only positions.

    :return: A list of MultipleSeqAlignment objects, one for each desired output
        size.
    """
    records = [r for r in alignment]
    sorted_dr_pairs = sorted(zip(records, distances),
                             key=operator.itemgetter(1))
    sorted_records = map(operator.itemgetter(0), sorted_dr_pairs)
    record_sets = []
    for size in sizes:
        rset = sorted_records[:size]
        record_sets.append(rset)

    if contract is True:
        alns = [contract_alignment(MultipleSeqAlignment(rset))
                for rset in record_sets]
    else:
        alns = [MultipleSeqAlignment(rset) for rset in record_sets]

    return alns


def mix_closest_and_furthest(alignment, distances, k_closest, k_furthest,
                             contract=True):
    """Returns an alignment from the ``k_closest`` and ``k_furthest`` records
    from the input alignment based on the given distances.

    :param alignment:
    :param distances:
    :param k_closest:
    :param k_furthest:
    :param contract:
    :return:
    """
    if (k_closest + k_furthest) > len(alignment):
        raise ValueError('Too many sequences to choose, would overlap.'
                         ' (closest: {0}, furhtest: {1}, available: {2})'
                         ''.format(k_closest, k_furthest, len(alignment)))
    records = [r for r in alignment]
    sorted_dr_pairs = sorted(zip(records, distances),
                             key=operator.itemgetter(1))
    sorted_records = map(operator.itemgetter(0), sorted_dr_pairs)
    closest_records = sorted_records[:k_closest]
    furthest_records = sorted_records[-k_furthest:]
    output = MultipleSeqAlignment(itertools.chain(closest_records,
                                                  furthest_records))
    return output

#
# def plot_distances(alignment, pivot, bins=50):
#     """Plots the distances of records in the given alignment with respect to the
#     pivot sequence. Plots both the individual distances and a histogram using
#     the given number of bins (50 by default).
#     """
#     distances = distances_to_seq(alignment, pivot)
#     plt.figure(1)
#     plt.subplot(211)
#     plt.plot(distances)
#     plt.subplot(212)
#     plt.hist(distances, bins)
#     plt.show()


def generate_distributed_aln_names(bins, pivot_name, tprof_name):
    """Generates names for distributed alignments according to the given bins.
    The naming scheme assumes four types of distributed alignments: cumulative
    and noncumulative, and both with and without pivot.

    * With the pivot sequence:
        * All       ``pivot.i.tprofname.aln``
        * By range of distances from pivot sequence: ``pivot.i.tprofname.min-max.aln``
        * By range, cumulative: ``pivot.i.tprofname.-max.aln``
    * Without pivot sequence:
        * By range of distances from pivot sequence: ``pivot.tprofname.min-max.aln``
        * By range, cumulative ``pivot.tprofname.-max.aln``

    :param bins:

    :param pivot_name:

    :param tprof_name:

    :return: cnames_nonpivot, cnames_withpivot, names_nonpivot, names_withpivot
    """
    logging.info(
        'Generating alignment names from bins {0}...'.format(bins))
    names_nonpivot = ['.'.join([pivot_name,
                                tprof_name,
                                '{0}-{1}'.format(bins[i],
                                                 bins[i + 1]),
                                'aln']) for i, _ in enumerate(bins[:-1])]
    names_withpivot = ['.'.join([pivot_name,
                                 'w',
                                 tprof_name,
                                 '{0}-{1}'.format(bins[i],
                                                  bins[i + 1]),
                                 'aln']) for i, _ in enumerate(bins[:-1])]
    cnames_nonpivot = ['.'.join([pivot_name,
                                 tprof_name,
                                 '{0}-{1}'.format(bins[0], bins[i + 1]),
                                 'aln']) for i, _ in enumerate(bins[:-1])]
    cnames_withpivot = ['.'.join([pivot_name,
                                  'w',
                                  tprof_name,
                                  '{0}-{1}'.format(bins[0], bins[i + 1]),
                                  'aln']) for i, _ in enumerate(bins[:-1])]
    return cnames_nonpivot, cnames_withpivot, names_nonpivot, names_withpivot


def generate_closest_aln_names(sizes, pivot_name, tprof_name):
    """Generates the names for fixed-size closest-sequence alignments. There
    are two kinds: with the pivot and without the pivot, and the alignments
    are always cumulative.

    >>> generate_closest_aln_names([1, 10], 'p', 'tprof')
    (['p.tprof.closest_1.aln', 'p.tprof.closest_10.aln'], ['p.w.tprof.closest_1.aln', 'p.w.tprof.closest_10.aln'])
    """
    logging.info('Generating alignment names from sizes {0}...'.format(sizes))
    names_nonpivot = ['.'.join([pivot_name,
                                tprof_name,
                                'closest-{0}'.format(s),
                                'aln']) for s in sizes]
    names_withpivot = ['.'.join([pivot_name,
                                 'w',
                                 tprof_name,
                                 'closest-{0}'.format(s),
                                 'aln']) for s in sizes]
    return names_nonpivot, names_withpivot


def generate_closest_and_furthest_aln_names(k_closest, k_furthest, pivot_name,
                                            tprof_name):
    logging.info('Generating alignment names from closest {0}, furthest {1}'
                 '...'.format(k_closest, k_furthest))
    names_nonpivot = ['.'.join([pivot_name,
                                tprof_name,
                                'c{0}-f{1}'.format(c, f),
                                'aln']) for c, f in zip(k_closest, k_furthest)]
    names_withpivot = ['.'.join([pivot_name,
                                 'w',
                                 tprof_name,
                                 'c{0}-f{1}'.format(c, f),
                                 'aln']) for f, c in zip(k_furthest, k_closest)]
    return names_nonpivot, names_withpivot


def restore_description_by_id(records_1, records_2):
    """For each sequence record in records_1, will find a record from record_2
    with the same ID and output the list of these records from records_2.

    This is a utility function to help extract taxonomy information in record
    ``description`` members but gets lost when the set of records gets aligned.
    """
    records_by_id = {rec.id: rec for rec in records_2}
    for rec in records_1:
        rec.description = records_by_id[rec.id].description


def species_from_alignment(alignment_file, seqfile_with_taxonomy):
    """When aligning sequences in Biopython through the ClustalwCommandLine,
    the ``description`` field of the sequences gets lost. This field contains
    taxonomy information that we'd have liked to retain. So, we take the aligned
    records and paste the description back into them based on the original
    sequence file with the taxonomy information."""
    aln = read_alignment(alignment_file)
    seqs = [seq for seq in SeqIO.parse(open(seqfile_with_taxonomy), 'fasta')]
    aln_records = [copy.deepcopy(r) for r in aln]
    restore_description_by_id(aln_records, seqs)
    output_alignment = MultipleSeqAlignment(aln_records)
    return output_alignment