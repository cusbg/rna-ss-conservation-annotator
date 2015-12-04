#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains constants for the PseudoknotSecstruc class operation.

The constants are symbols used in dot-paren strings with pseudoknots.
"""
__author__  = "Jan Hajic"
__contributors__ = ""
__credits__ = [""]
__maintainer__ = "Jan Hajic"
__email__ = "hajicj@ufal.mff.cuni.cz"
__status__ = "Prototype"

#: a dictionary of matching brackets for base pairs in dot-paren notation:
#: ``( => )``, ``[ => ]``, etc.
symbols = {
    '(':')', '[':']', '{':'}', '<':'>', 
    'a':'A', 'b':'B', 'c':'C', 'd':'D', 'e':'E', 'f':'F', 
    'g':'G', 'h':'H', 'i':'I', 'j':'J', 'k':'K', 'l':'L', 
    'm':'M', 'n':'N', 'o':'O', 'p':'P', 'q':'Q', 'r':'R', 
    's':'S', 't':'T', 'u':'U', 'v':'V', 'w':'W', 'x':'X', 'y':'Y', 'z':'Z', 
    }

#: A list of base pair start symbols: ``(``, ``[``, etc.
open_symbols = symbols.keys()

#: A list of base pair end symbols: ``)``, ``]``, etc.
close_symbols = symbols.values()
    
#: A dictionary of matching brackets in dot-paren notation that are used
#: to indicate pseudoknotted base pairs.
pknot_symbols = { o : symbols[o] for o in symbols if o != '(' }

#: A list of pseudoknotted base pair start symbols: ``[``, ``{``, etc.
pknot_open_symbols = pknot_symbols.keys()

#: A list of pseudoknotted base pair end symbols: ``]``, ``}``, etc.
pknot_close_symbols = pknot_symbols.values()

#: A set of characters that are not pseudoknotted: ``(``, ``)`` and ``.``
non_pknot_chars = {'(', ')', '.'}

#: ``bracket_priorities``: A set of opening symbols in the order of
#: pseudoknottedness level (first member is ``(``, then ``[``, etc.)
bracket_priorities = ['(', '[', '{', '<', 'a', 'b', 'c', 'd', 'e']
