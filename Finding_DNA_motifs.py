"""
Problem (see here: https://rosalind.info/problems/subs/)
Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s
(as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i
of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s
; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s
(see the Sample below).

Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
"""

import re

sequence_1 = 'GATATATGCATATACTT'
sequence_2 = 'ATAT'

def motif_function(my_string, my_substring):
    assert len(my_substring) <= len(my_string)

    matches = re.finditer(my_substring, my_string)

    list_of_match_start_indices = []

    for match in matches:
        start = match.start() + 1 # adding one to translate between python indexing and real world indexing
        list_of_match_start_indices.append(start)

    return list_of_match_start_indices

