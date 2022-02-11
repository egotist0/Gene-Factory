from random import *
import numpy as np

pairwise = "AAATTTCCGG"


def max(a, b):
    if a >= b:
        return a
    else:
        return b


def Pairwise_alignment(len1, len2, str1, str2):
    Score = np.zeros((200, 200))
    gap = -2.5
    m = 5
    mm = -4

    for i in range(len1 + 1):
        Score[i, 0] = i * gap
    for j in range(len2 + 1):
        Score[0, j] = j * gap

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if str1[i - 1] == str2[j - 1]:
                Score[i, j] = Score[i - 1, j - 1] + m
            else:
                Score[i, j] = Score[i - 1, i - 1] + mm
            if Score[i, j] < max(Score[i - 1, j], Score[i, j - 1]) + gap:
                Score[i][j] = max(Score[i - 1, j], Score[i, j - 1]) + gap
    return Score


def printAlign(sc, i, j, s1, s2, saln, raln):
    if not (i and j):
        return ''.join(saln), ''.join(raln)
    if sc[i - 1, j] >= sc[i, j - 1] and sc[i - 1, j] >= sc[i - 1, j - 1]:
        saln.append(s1[i - 1])
        raln.append('-')
        printAlign(sc, i - 1, j, s1, s2, saln, raln)
    elif sc[i - 1, j - 1] >= sc[i, j - 1] and sc[i - 1, j - 1] >= sc[i - 1, j]:
        saln.append(s1[i - 1])
        raln.append(s2[j - 1])
        printAlign(sc, i - 1, j - 1, s1, s2, saln, raln)

    else:
        saln.append('-')
        raln.append(s2[j - 1])
        printAlign(sc, i, j - 1, s1, s2, saln, raln)


def main(st1, st2):
    len1 = len(st1)
    len2 = len(st2)
    score = Pairwise_alignment(len1, len2, st1, st2)
    l1 = []
    l2 = []
    printAlign(score, len1, len2, st1, st2, l1, l2)
    return score[len1, len2], ''.join(l1[::-1]), ''.join(l2[::-1])
