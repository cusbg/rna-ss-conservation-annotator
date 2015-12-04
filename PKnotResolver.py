#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# PKnotResolver.py
"""
This module handles selecting which of conflicting base pairs are the
pseudoknots and which belong to the "un-knotted" part.

Selecting pseudoknots is a non-trivial task. For instance, which of the
following solutions should be used?::

  ..(((..[[)))......]]
  ..[[[..((]]]......]]

We would choose the first one if we want to minimize the number of pseudoknotted
base pairs. We would choose the second one if we want to minimize pseudoknot
length.

The algorithm in this module is a greedy graph coloring algorithm
using a Breadth-First Search with a specific initialization and coloring
priority strategy.

-----------------------

"""
__author__ = "Jan Hajic"
__contributors__ = ""
__credits__ = ["Michal Klimpera", "Jan Drozen"]
__maintainer__ = "Jan Hajic"
__email__ = "hajicj@ufal.mff.cuni.cz"
__status__ = "Prototype"

import logging


class UndyingIndexedQueue(object):
    """A fast, not very memory-efficient queue for the BFS search of the
    helix conflict graph. Doesn't delete on pop; only moves an index to the
    next position in the list of elements push()-ed so far.

    The queue is *indexed*, meaning that it is possible to quickly check
    whether a given item is in the queue or not.

    Undefined behavior occurrs if an item occurrs more than once into the
    "active" part of the queue: after popping the first occurrence, the index
    will delete the item, but the item will still be in the queue somewhere
    later on. Before push()-ing, always check by has(). Since this queue is
    intended for a BFS graph search, this restriction is sensible: no vertex
    should be pushed into the queue more than once per search.

    Helper class for :class:`PKnotResolver`.
    """
    def __init__(self):
        self.q = []
        self._index = set(self.q)
        self.out = 0

    def push(self, p):
        """Adds an item to the end of the queue.

        :type p: Any
        :param p: The item to add.
        """
        if p in self._index:
            logging.warn('Adding already existing item to unique queue! Undefined behavior will result.')
        self.q.append(p)
        self._index.add(p)

    def pop(self):
        """Retrieves an item from the start of the queue.

        :rtype: Any
        :returns: The item at the start of the queue.
        """
        if self.out >= len(self.q):
            return None
        k = self.q[self.out]
        if k not in self.q[self.out+1:]:
            self._index.remove(k)
        self.out += 1
        return k

    def has(self, p):
        """Checks whether ``p`` is already in the queue.

        :type p: Any
        :param p: An item which is checked for existence in the queue.

        :rtype: Boolean
        :returns: True on item ``p`` found, False if not found.
        """
        return p in self._index

    def isempty(self):
        return len(self.q) == 0 or self.out >= len(self.q)


class PKnotResolver(object):
    """This class resolves which base pairs should be called
    pseudoknots and which shouldn't.

    Assigns levels of "pseudoknottedness" using a greedy BFS graph-coloring
    algorithm.
    """
    
    colors = ['(', '[', '{', '<']

    def __init__(self):
        """This a static class - no initiation necessary, all instances will
        be the same."""
        pass

    def color_bplist(self, base_pairs):
        """Given a list of base pairs, returns a list of base pairs
        for each level of pseudoknottedness."""

        if base_pairs is None or base_pairs == []:
            return {}

        # Algorithm:
        helices, bp2helix = self.bplist2helixlist(base_pairs)
        # helices is just a minimized base pair list: each helix is
        # represented by one base pair. Each helix also has some properties:
        # length, number of base pairs, number of spanned helices, etc.

        #  - build helix conflict graph
        hc_graph = self.build_helix_conflict_graph(helices)
        # The hc_graph is a graph in adjacency list format (list of lists).

        #logging.debug('Helices:\n%s' % '\n'.join(str(h['bps'][0]) for h in helices))
        #logging.debug('Helix conflict graph:\n%s' % '\n'.join([str(h) for h in hc_graph ]))

        #  - color helix conflict graph
        helix_coloring, helix_coloring_index = self.color_bipartite_hc_graph(hc_graph, helices)
        # helix_coloring is a dict of helix IDs, with 'color' keys
        # ('(', '.', etc.) and helix ID values

        bp_coloring = { color : self.helixlist2bplist(helix_coloring[color], helices) for color in helix_coloring }

        return bp_coloring

    def bplist2helixlist(self, bplist):
        """Given a list of base pairs, returns a list of representant
        base pairs for each helix (one per helix). Each member of the list
        is a dictionary with additional information about the helix (number of
        base pairs, number of spanned helices, etc.).

        A helix is defined as a set of two regions where the i-th residue
        from the 5'-end of one is paired to the i-th residue from the 3'-end
        of the other. This is a helix::

          ((()))

        This is not a helix, but two helices::

          ((.()))

        :type bplist: list(tuple(int, int))
        :param bplist: A list of base pairs. Each base pair is represented as
            as a tuple of indices of participating residues (indices in some
            source string, not PDB indices).

        :rtype: list(dict( 'bps' => list(tuple(int, int)), 'size' => int,
            'length' => int, 'spanning' => int'))
        :returns: A list of Helix dict items. Each Helix item has the following
            members:

            * ``bps`` : a list of base pairs participating in the helix

            * ``size`` : the number of base pairs participating in the helix

            * ``length`` : the length of the longest (outermost) base pair in
              the helix

            * ``spanning`` : the number of helices nested within this helix
              [NOT IMPLEMENTED]
        """
        if not bplist:
            return []

        if len(bplist) == 1:
            return [ {'bps' : [bplist[0]], 
                      'size' : 1, 
                      'length' : bplist[0][1] - bplist[0][0] + 1, 
                      'spanning' : 1 } ], { bplist[0] : 0 }

        # Just in case the list is not sorted on input.
        sorted_bplist = sorted(bplist)

        helices = []

        prev_bp = sorted_bplist[0]
        current_helix_bps = [prev_bp]
        for bp in sorted_bplist[1:]:

            if bp[0] == prev_bp[0] + 1 and bp[1] == prev_bp[1] - 1:
                current_helix_bps.append(bp)
            else:
                helix = { 'bps' : current_helix_bps, 
                          'size' : len(current_helix_bps),
                          'length' : current_helix_bps[0][1] - current_helix_bps[0][0] + 1,
                          'spanning' : 1 } # Can't really compute spanning yet.
                helices.append(helix)

                current_helix_bps = [bp]

            prev_bp = bp

        # Last base pair
        if current_helix_bps:
            helix = { 'bps' : current_helix_bps, 
                      'size' : len(current_helix_bps),
                      'length' : current_helix_bps[0][1] - current_helix_bps[0][0] + 1,
                      'spanning' : 1 } # Can't really compute spanning yet.
            helices.append(helix)

        bp2helix = {} # For each base pair (tuple), the index of the helix it belongs to
        for h_idx, helix in enumerate(helices):
            for bp in helix['bps']:
                bp2helix[bp] = h_idx            

        return helices, bp2helix

    def helixlist2bplist(self, h_idx_list, helices):
        """From a given list of helix indices and the helices list,
        reconstructs the sorted list of base pairs belonging to the indexed
        helices.

        :type h_idx_list: list(int)
        :param h_idx_list: A list of indices into the ``helices`` parameter.

        :type helices: list(helix)
        :param helices: A list of helix dicts (see ``bplist2helixlist``).

        :rtype: list(tuple(int, int))
        :returns: A list of base pairs that participate in the helices given
            by ``h_idx_list``. The base pairs are obtained from the ``bps``
            members of the helix dicts.
        """
        bplist = []
        for h_idx in h_idx_list:
            bplist.extend(helices[h_idx]['bps'])
        return sorted(bplist)

    def _is_conflict(self, h, h_prime):
        """Are the given two helices interlocking?

        :type h: helix dict
        :param h: The first helix.

        :type h_prime: helix dict
        :param h_prime: The second helix.

        :rtype: Boolean
        :returns: True if the given helices are in conflict (i.e. are not
            well-nested). False otherwise.
        """
        bp = h['bps'][0]
        bp_prime = h_prime['bps'][0]

        return (bp[0] < bp_prime[0] and bp[1] > bp_prime[0] and bp_prime[1] > bp[1]) or (bp_prime[0] < bp[0] and bp_prime[1] > bp[0] and bp[1] > bp_prime[1])

    def build_helix_conflict_graph(self, helices):
        """Given a list of helices, builds a graph in adjacency list
        representation. A helix is a (start, end) tuple, like a base pair;
        the graph indexes them in the order they are given.

        :type helices: list(helix dict)
        :param helices: A list of helix dicts.

        :rtype: list(list(int))
        :returns: A conflict graph in adjacency list format. The i-th member of
            the top-level list represents the list of helices conflicting with
            the i-th helix from the ``helices`` parameter.

        """
        hc_graph = [ [] for h in helices ]

        for i, h in enumerate(helices):

            for j, h_prime in enumerate(helices[i+1:]):

                if self._is_conflict(h, h_prime):
                    hc_graph[i].append(i + j + 1)
                    hc_graph[i+j+1].append(i)

        return hc_graph

    def color_bipartite_hc_graph(self, hc_graph, helices):
        """Given a graph in adjacency list format, attempts to color it
        greedily with an ordering specified by the sort_helices() method.

        The "colors" used in the coloring are different bracket types.

        This method implements the core algorithm of PKnotResolver.

        :type hc_graph: list(list(int))
        :param hc_graph: The conflict graph of the ``helices``, as computed
            by the ``build_helix_conflict_graph`` method.

        :type helices: list(helix dict)
        :param helices: The list of helices, as computed by the
            ``bplist2helixlist`` method from the original base pair list.

        :rtype: tuple(dict(color => list(helix_idx)), dict(helix_idx => color))
        :returns: The ``coloring`` dict and the ``coloring_index`` dict.

            The ``coloring`` dict is a dictionary that for each color assigned
            to at least one helix returns the list of helices (indices into the
            ``helices`` parameter) colored with this particular color.

            The ``coloring_index`` dict stores for each helix index the color
            with which the given helix was colored.
        """

        if self._no_graph_conflicts(hc_graph):
            coloring = { '(': range(len(helices)) }
            coloring_idx = { h: '(' for h in xrange(len(hc_graph)) }
            return coloring, coloring_idx

        ordering = self.order_helices(hc_graph, helices)

        #logging.debug('Ordering: %s' % str(ordering))

        coloring = { '(': [], '[': [] }
        coloring_index = {}
        color_priorities = [ '(', '[', '{', '<', 'a', 'b', 'c', 'd', 'e' ]
              # There shouldn't be more colors than this.

        uncolored = set(ordering)
        colored = set([])

        queue = UndyingIndexedQueue()
        #logging.debug('Queue construction state: %s, with out at %d' % (str(queue.q), queue.out))

        queue.push(ordering[0])

        #logging.debug('Queue initial state: %s, with out at %d' % (str(queue.q), queue.out))

        is_first_in_component = True
        is_first_in_nontrivial_component = False
        while not queue.isempty():

            #logging.debug('Current queue state: %s' % str(queue.q[queue.out:]))
            #logging.debug('Current coloring state:\n%s' % str(coloring))

            h_idx = queue.pop()

            #logging.debug('Current h_idx: %d' % h_idx)

            neighbors = hc_graph[h_idx]

            is_first_in_nontrivial_component = is_first_in_component and neighbors

            # Find color to use
            if is_first_in_nontrivial_component:
                current_color = color_priorities[1]

                #logging.debug('Coloring first of nontrivial component: %d' % h_idx)

                is_first_in_component = False
            else:
                color_iter = color_priorities.__iter__()
                current_color = color_iter.next()
                for neighbor in neighbors:
                    if neighbor in colored:
                        if current_color == coloring_index[neighbor]:
                            current_color = color_iter.next()
                            
                            if current_color in color_priorities[2:]:
                                logging.warn('The conflict graph is NOT bipartite! At: helix %d, represented as %s' % (h_idx, str(helices[h_idx]['bps'][0])))

                #logging.debug('Coloring regular: %d' % h_idx)

            # Color vertex
            if current_color not in coloring:
                coloring[current_color] = [h_idx]
            else:
                coloring[current_color].append(h_idx)

            coloring_index[h_idx] = current_color

            # Update coloration tracking
            colored.add(h_idx)
            uncolored.remove(h_idx)

            # Append to stack - in given order
            ordered_neighbors =  self.order_helices([hc_graph[h] for h in neighbors], [ helices[h] for h in neighbors], neighbors)

            #logging.debug('Ordered neighbors: %s' % str(ordered_neighbors))

            for neighbor in ordered_neighbors:
                if neighbor in uncolored and not queue.has(neighbor):
                    queue.push(neighbor)

            # If we reached an end of a connected component:
            # choose next in ordering
            if queue.isempty():
                is_first_in_component = True
                o_iter = ordering.__iter__()
                h_idx_candidate = o_iter.next()
                try:
                    while h_idx_candidate in colored:
                        h_idx_candidate = o_iter.next()
                except StopIteration:
                    continue
                queue.push(h_idx_candidate)

        return coloring, coloring_index

    def order_helices(self, hc_graph, helices, hc_graph_indices = None):
        """Returns a list of indices into the list of helices that is ordered
        according to some criteria. In the case of pseudoknot assignment,
        this criteria is the vertex degree (number of conflicting helices).

        :type hc_graph: list(list(int))
        :param hc_graph: The conflict graph from which the indices are taken.

        :type helices: list(helix dict)
        :param helices: The list of helices from which the conflict graph is
            generated.

        :type hc_graph_indices: list(int)
        :param hc_graph_indices: A list of indices into both ``hc_graph`` and
            ``helices`` that denotes which conflict graph vertices (= helices)
            should be ordered.

        :rtype: list(int)
        :returns: A list of the ``hc_graph_indices`` sorted in order of
            priority: vertices that come first in the list will be colored
            first.
        """

        if hc_graph_indices is None:
            hc_graph_indices = range(len(hc_graph))

        return [ hc_graph_indices[hc_g_idx] for hc_g_idx in sorted(range(len(hc_graph)), key=lambda k: len(hc_graph[k]), reverse=True) ]

    def _no_graph_conflicts(self, hc_graph):
        """Checks whether there is at least one edge in the conflict graph.

        :type hc_graph: list(list(int))
        :param hc_graph: The conflict graph to check for edge existence.

        :rtype: Boolean
        :returns: True if at least one conflict is found (at least one edge in
            conflict graph - one non-empty neighbor list)

        """
        no_conflict = True
        for h_edges in hc_graph:
            if h_edges:
                no_conflict = False
                break

        return no_conflict
