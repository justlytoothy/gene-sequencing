#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the 3 for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class Node:
    val = 0
    prev = None

    def __init__(self, val, prev):
        self.val = val
        self.prev = prev


class GeneSequencing:
    alignedSeq1 = ''
    alignedSeq2 = ''
    optimalScore = 0

    def __init__(self):
        self.alignedSeq1 = ''
        self.alignedSeq2 = ''
        self.optimalScore = 0
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        ###################################################################################################
        # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
        if banded:
            self.banded_alignment(seq1, seq2)
        else:
            self.unbanded_alignment(seq1, seq2)

        score = self.optimalScore
        # TODO add aligned seqs here
        alignment1 = '{}  DEBUG:({} chars,align_len={}{})'.format(
            self.alignedSeq1, len(seq1), align_length, ',BANDED' if banded else '')
        alignment2 = '{}  DEBUG:({} chars,align_len={}{})'.format(
            self.alignedSeq2, len(seq2), align_length, ',BANDED' if banded else '')
        ###################################################################################################

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}

    # O(kn) max time complexity because of the O(kn) loop to fill
    # and other O(kn) loop to reconstruct alignment strings
    # Max space of kn because of the map that gets filled with k * n items
    def banded_alignment(self, seq1, seq2):
        # Initialize the map
        dp_map = {}
        # Chop off anything beyond max align length
        len_one = min(self.MaxCharactersToAlign, len(seq1))
        len_two = min(self.MaxCharactersToAlign, len(seq2))
        self.alignedSeq1 = ""
        self.alignedSeq2 = ""
        # Check if alignment is possible
        if abs(len_one - len_two) > MAXINDELS:
            self.alignedSeq1 = "No Alignment Possible."
            self.alignedSeq1 = "No Alignment Possible."
            self.optimalScore = float('inf')
            return
        # Constant time because will only loop for the number of MAXINDELS + 1 which is constant
        for i in range(0, MAXINDELS + 1):
            if (i - 1, 0) in dp_map:
                coord = (i - 1, 0)
            else:
                coord = (-1, -1)
            dp_map[(i, 0)] = (coord[0], coord[1], i * INDEL)
        # Constant time because will only loop for the number of MAXINDELS + 1 which is constant
        for i in range(0, MAXINDELS + 1):
            if (0, i - 1) in dp_map:
                coord = (0, i - 1)
            else:
                coord = (-1, -1)
            dp_map[(0, i)] = (coord[0], coord[1], i * INDEL)
        # O(kn) because loops n times with inner loop that loops k times
        # This fills the dp_map with k * n items, giving it space of kn
        for i in range(1, len_one + 1):             # O(n) because loops for len_one
            x = i - 1
            for j in range(i - MAXINDELS, i + 1 + MAXINDELS):           # O(k) where k is MAXINDELS * 2 + 1
                y = j - 1
                one, two, three = float('inf'), float('inf'), float('inf')
                # All accessing of the dp_map is constant because of the data structure used
                if 0 < j < len_two + 1:
                    if (x, y) in dp_map:
                        one = (MATCH if seq1[x] == seq2[y] else SUB) + dp_map[(x, y)][2]
                    if (i, y) in dp_map:
                        two = INDEL + dp_map[(i, y)][2]
                    if (x, j) in dp_map:
                        three = INDEL + dp_map[(x, j)][2]
                    # min operation is of inconsequential time
                    min_val = min(one, two, three)
                    if two == min_val:
                        coord = (i, y)
                    elif three == min_val:
                        coord = (x, j)
                    else:
                        coord = (x, y)
                    # Constant time to assign
                    dp_map[(i, j)] = (coord[0], coord[1], min_val)
        curr = (len_one, len_two)
        prev = (dp_map[curr][0], dp_map[curr][1])
        # O(kn) time maximum, kn is size of the dp_map
        while prev != (-1, -1):
            if prev == (curr[0], curr[1] - 1):
                self.alignedSeq1 = "-" + self.alignedSeq1
                self.alignedSeq2 = seq2[curr[1] - 1] + self.alignedSeq2
            elif prev == (curr[0] - 1, curr[1]):
                self.alignedSeq1 = seq1[curr[0] - 1] + self.alignedSeq1
                self.alignedSeq2 = "-" + self.alignedSeq2
            elif prev == (curr[0] - 1, curr[1] - 1):
                self.alignedSeq1 = seq1[curr[0] - 1] + self.alignedSeq1
                self.alignedSeq2 = seq2[curr[1] - 1] + self.alignedSeq2
            curr = prev
            prev = (dp_map[curr][0], dp_map[curr][1])
        # All constant time assignments
        self.alignedSeq1 = self.alignedSeq1[:100]           # Chopping to first 100
        self.alignedSeq2 = self.alignedSeq2[:100]
        self.optimalScore = dp_map[(len_one, len_two)][2]

    # O(mn) max time complexity
    # Max space of mn because of the map and array that gets filled with m * n items
    def unbanded_alignment(self, seq1, seq2):
        # i is word on left, j is word on top
        len_one = min(self.MaxCharactersToAlign, len(seq1))
        len_two = min(self.MaxCharactersToAlign, len(seq2))
        self.alignedSeq1 = ""
        self.alignedSeq2 = ""
        back_map = {}
        # Time of O(mn) also where arr gets filled with m * n items
        arr = [[0 for i in range(len_two + 1)] for j in range(len_one + 1)]
        # O(n)
        for i in range(len_one + 1):
            arr[i][0] = i * INDEL
            back_map[(i, 0)] = 'top'
        # O(m)
        for i in range(len_two + 1):
            arr[0][i] = i * INDEL
            back_map[(0, i)] = 'left'
        # O(mn) to loop n times with loop for m times inside
        # Also where arr and map both get filled again to max of m * n space
        for i in range(1, len_one + 1):                 # O(n)
            x = i - 1
            for j in range(1, len_two + 1):             # O(m)
                y = j - 1
                one = (MATCH if seq1[x] == seq2[y] else SUB) + arr[x][y]
                two = INDEL + arr[i][y]
                three = INDEL + arr[x][j]
                min_val = min(one, two, three)
                arr[i][j] = min_val
                if two == min_val:
                    back_map[(i, j)] = 'left'
                elif three == min_val:
                    back_map[(i, j)] = 'top'
                else:
                    back_map[(i, j)] = 'diagonal'

        i, j = len_one, len_two
        # O(mn) maximum time here as well
        while i > 0 or j > 0:
            if back_map[(i, j)] == 'diagonal':
                self.alignedSeq1 = seq1[i - 1] + self.alignedSeq1
                self.alignedSeq2 = seq2[j - 1] + self.alignedSeq2
                i -= 1
                j -= 1
            elif back_map[(i, j)] == 'top':
                self.alignedSeq1 = seq1[i - 1] + self.alignedSeq1
                self.alignedSeq2 = '-' + self.alignedSeq2
                i -= 1
            else:
                self.alignedSeq1 = '-' + self.alignedSeq1
                self.alignedSeq2 = seq2[j - 1] + self.alignedSeq2
                j -= 1
        self.alignedSeq1 = self.alignedSeq1[:100]
        self.alignedSeq2 = self.alignedSeq2[:100]
        self.optimalScore = arr[len_one][len_two]
