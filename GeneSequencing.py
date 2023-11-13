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

# Used to compute the bandwidth for banded version
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

    def banded_alignment(self, seq1, seq2):
        pass

    def unbanded_alignment(self, seq1, seq2):
        # i is word on left, j is word on top
        len_one = min(self.MaxCharactersToAlign, len(seq1))
        len_two = min(self.MaxCharactersToAlign, len(seq2))
        self.alignedSeq1 = ""
        self.alignedSeq2 = ""
        arr = [[(0, 0, 0) for i in range(len_two + 1)] for j in range(len_one + 1)]
        curr = 0

        prev = 0
        for i in range(len_one + 1):
            arr[i][0] = (prev, 0, curr)
            prev = i
            curr += 5
        curr = 0
        prev = 0
        for i in range(len_two + 1):
            arr[0][i] = (0, prev, curr)
            prev = i
            curr += 5
        for i in range(1, len_one + 1):
            x = i - 1
            for j in range(1, len_two + 1):
                y = j - 1
                one = (-3 if seq1[x] == seq2[y] else 1) + arr[i - 1][j - 1][2]
                two = 5 + arr[i][j - 1][2]
                three = 5 + arr[i - 1][j][2]
                min_val = min(one, two, three)
                if one == min_val:
                    val = (i - 1, j - 1, min_val)
                elif two == min_val:
                    val = (i, j - 1, min_val)
                else:
                    val = (i - 1, j, min_val)
                arr[i][j] = val
        curr = arr[len_one][len_two]
        i, j = len_one - 1, len_two - 1
        while i >= 0 or j >= 0:
            next_val = arr[curr[0]][curr[1]]
            same_x = next_val[0] == curr[0]
            same_y = next_val[1] == curr[1]
            if same_x and same_y:
                self.alignedSeq1 = seq1[i] + self.alignedSeq1
                self.alignedSeq2 = seq2[j] + self.alignedSeq2
                # if i >= 0:
                #     self.alignedSeq1 = seq1[i] + self.alignedSeq1
                # else:
                #     self.alignedSeq1 = '-' + self.alignedSeq1
                # if j >= 0:
                #     self.alignedSeq2 = seq2[j] + self.alignedSeq2
                # else:
                #     self.alignedSeq2 = '-' + self.alignedSeq2
                i -= 1
                j -= 1
            elif same_y:
                self.alignedSeq1 = seq1[i] + self.alignedSeq1
                # if i >= 0:
                #     self.alignedSeq1 = seq1[i] + self.alignedSeq1
                # else:
                #     self.alignedSeq1 = '-' + self.alignedSeq1
                self.alignedSeq2 = '-' + self.alignedSeq2
                i -= 1
            elif same_x:
                self.alignedSeq1 = '-' + self.alignedSeq1
                self.alignedSeq2 = seq2[j] + self.alignedSeq2
                j -= 1
            curr = next_val
        self.optimalScore = arr[len_one][len_two][2]
