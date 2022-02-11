from random import *
import numpy as np
import matplotlib.pyplot as pl
from scipy import stats

pairwise = "AAATTTCCGG"
score = np.zeros((200, 200))
SCORE = open("score.txt", 'w')
SEQ = open("seq.txt", 'w')
OUT = open("outcome.txt",'w')

def getseq():
    seq = []
    num = randint(40, 60)
    for i in range(num):
        seq.append(pairwise[randint(0, 9)])
    return seq


def Pairwise_alignment(str1, str2):
    len1 = len(str1)
    len2 = len(str2)
    score = np.zeros((200, 200))
    st1 = []
    st2 = []
    gap = -2.5
    m = 5
    mm = -4

    for i in range(len1+1):
        score[i, 0] = i * gap
    for j in range(len2+1):
        score[0, j] = j * gap

    for i in range(1, len1+1):
        for j in range(1, len2+1):
            if str1[i - 1] == str2[j - 1]:
                score[i, j] = score[i - 1, j - 1] + m
            else:
                score[i, j] = score[i - 1, i - 1] + mm
            if score[i, j] < max(score[i - 1, j], score[i, j - 1]) + gap:
                score[i][j] = max(score[i - 1, j], score[i, j - 1]) + gap
    SEQ.write(''.join(str1) + " 和\n" + ''.join(str2) + "\n的序列比对的得分为" + str(score[len1, len2]))
    SCORE.write(str(score[len1, len2]) + '\n')
    printAlign(score, len1, len2, str1, str2, st1, st2)
    SEQ.write("\n" * 2)


def printAlign(score, i, j, s1, s2, saln, raln):
    if not (i and j):
        SEQ.write("\n最佳匹配结果为\n" + ''.join(saln) + " 和\n" + ''.join(raln))
        return 0
    if score[i - 1, j] >= score[i, j - 1] and score[i - 1, j] >= score[i - 1, j - 1]:
        saln.append(s1[i - 1])
        raln.append('-')
        printAlign(score, i - 1, j, s1, s2, saln, raln)
    elif score[i - 1, j - 1] >= score[i, j - 1] and score[i - 1, j - 1] >= score[i - 1, j]:
        saln.append(s1[i - 1])
        raln.append(s2[j - 1])
        printAlign(score, i - 1, j - 1, s1, s2, saln, raln)

    else:
        saln.append('-')
        raln.append(s2[j - 1])
        printAlign(score, i, j - 1, s1, s2, saln, raln)


def loadData():
    fp = open("score.txt", 'r')
    lines = fp.readlines()
    x = []

    for line in lines:
        line = line.replace("\n", "")
        line = float(line)
        x.append(line)
    return np.array(x)


def draw_hist(lenths):
    data = lenths
    # 对数据进行切片
    bins = np.linspace(min(data), max(data), 100)

    pl.hist(data, bins)
    pl.xlabel('Number of score')
    pl.ylabel('Number of occurences')
    pl.title('Frequency distribution of number of score')

    pl.savefig('outcome.png')
    pl.show()


def main():
    seq = []
    for i in range(50):
        seq.append(getseq())

    for i in range(50):
        for j in range(i, 50):
            Pairwise_alignment(seq[i], seq[j])
    SEQ.close()
    SCORE.close()
    data = loadData()
    draw_hist(data)
    pl.show()
    u = data.mean()  # 计算均值
    std = data.std()  # 计算标准差
    OUT.write("分数平均值为"+str(u)+'\n标准差为'+str(std)+'\nKS检验结果为：')
    OUT.write(str(stats.kstest(data, 'norm', (u, std))))
    OUT.close()
    SEQ.close()
    SCORE.close()
    return  stats.kstest(data, 'norm', (u, std))