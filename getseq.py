from random import *
pairwise = "AAATTTCCGG"
def getseq(num):
    seq = []
    for i in range(num):
        seq.append(pairwise[randint(0, 9)])
    return seq