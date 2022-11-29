# Gene comparison

Egotist





## I. Inspiration for the project

The inspiration for this project came from a semester-long course on *computational programming languages*, where the instructor was talking about dynamic programming algorithms that can be used for sequence matching to calculate the optimal score, which led to the question: **sequence matching**

> Sequence matching: is to determine the similarity or even homology between two or more sequences, and arrange them according to a certain law. Two or more sequences are aligned together and their similarities are marked. Spacers (usually indicated by a short horizontal line "-") can be inserted in the sequences.

This method is often used to study sequences that evolved from a common ancestor, especially biological sequences such as protein sequences or DNA sequences. In alignment, mismatches correspond to mutations, while null positions correspond to insertions or deletions. Sequence alignment can also be used for studies such as language evolution or inter-textual similarity.

Based on this problem, I adopted the Needleman/Wunsch algorithm (an algorithm dedicated to text matching) to solve the scoring problem of sequence matching and the best backtracking results, and then derived some small functions, such as random sequence generation and data analysis, etc. For ease of use consideration I designed a visual applet using the tinker package, and the whole program is in The whole program is carried out in anaconda3 environment of python3.8




## II. Main functions and brief introduction

Can be devided into 3 different components

### 1.Sequence matching and best traceback results

#### Ideas

In bioinformatics, we want to find some similarity relationship between two sequences, and this algorithm to find the similarity relationship of biological sequences is the double sequence comparison algorithm . The similarity between two sequences is usually determined by using the difference of characters between the two sequences. If the characters at the corresponding positions in the two sequences are different, then the similarity of the sequences is low, and vice versa, the similarity of the sequences is high.

- Biobase sequences consist of ATCG, a scoring rule needs to be set for comparison, I set `match=5,dismatch=-4,gap=-2.5` here, if the user needs to modify the rule, the location of the modified rule is in lines 26-28 of `data_analysis.py` file and lines 16-18 of `dp.py` file. Lines 16-18 of the `dp.py` file. ![](file:///home/egotist/PycharmProjects/dp_program/photo/1.jpg)

![](file:///home/egotist/PycharmProjects/dp_program/photo/2.jpg)



#### Scoring Functions

+ `gap` indicates a missing score of 2.5, `m` indicates a matching score of 5, and `mm` indicates a non-matching score set to -4.

+ Initialize the array, and for level 0, assign column 0 to the value `i*gap`;

+ The following two-layer `for` loop counts the `score` value for each position of the two-dimensional array.

  1. the overall per-site score is.

     **score for all three directions = score of the previous locus in that direction + score of the move**

     Finally, the highest score of the three directions is selected as the score of the bit, and the score of the whole matrix is obtained from top to bottom and from left to right in this cycle. 2.

  2. Then, after the previous step, we know that the score in the bottom right corner of ** must be the optimal score, because it is obtained from the optimal score of each seed case. **.

  3. The specific operation is

     1. if `str1[i - 1] = str2[j - 1]`, then `score[i, j]` is directly equal to `score[i-1, j-1]+m`
     2. if `str1[i - 1] ! = str2[j - 1]`, then `score[i, j]` is directly equal to `score[i-1, j-1]+mm`
     3. first carry out the above two steps, then the upper left corner `score[i, j-1]` is compared with the upper right corner `score[i-1, j]` to take out the larger value, then the larger value obtained is added to `gap`, if at this point `score[i, j]` is less than this new obtained score, then update `score[i, j]`, otherwise it is not updated.
     4. the final result is the selection of the highest score in the three directions as the score for that locus.

+ The score in the bottom right corner after the complete score table is formed is necessarily the best match score, which is written to the `SCORE` file; it is also returned to the applet

+ But we also need to know from which path it got this best score. So we need to backtrack, and here we call the `printAlign` function

#### 回溯函数

+ 回溯的方式就是看每个回溯位点的左上方，上方和左方最大值位置，最后就可以得到整个回溯路径；

+ `printAlign`函数输入的参数比较多，`score`是输入的二维评分数组，`i`和`j`表示的是此时回溯到的评分数组的位置，初始位置就是从(len1,len2)开始，`s1`和`s2`是进行的回溯的两条序列，`saln`和`raln`存储的是最后输出的最佳匹配结果

+ 三个`if`判断是是判断左上，左，上哪个得分最高，最高的就是最优的匹配路径：

  1. 从右下方开始，如果最大值出现在上面，则横向这条序列引入一个GAP ("-")，纵向这条序列取该处碱基；
  2. 如果最大值出现在左边，则纵向这条序列引入一个GAP ("-")，横向这条序列取该处碱基; 
  3. 如果最大值出现在左上角，则不引入GAP，纵向和横向均取该处碱基。
  4. 当然这个回溯的路径不一定唯一，当三个位置有两个相同的时候，两个路径都是可行的,但为了小程序输出方便，我用了`if...elif`这样就只会筛出一个序列输出

#### 使用手册

运行`smallprogram.py`文件，首先打开`On`，然后输入两条序列，点击上方按钮即可输出最高得分和最佳匹配方式

![](file:///home/egotist/PycharmProjects/dp_program/photo/%E5%88%9D%E5%A7%8B.jpg)

![](file:///home/egotist/PycharmProjects/dp_program/photo/3.png)

**由于这是一个优秀的文本字符串比对算法，因此不仅可以用来比对碱基序列，类似的文本比对，相似度比对，文字查重，蛋白质序列比对都可以同样进行**



### 2.生成比例为人类染色体的随机序列

#### 思路

人类1号染色体拥有2.3亿个碱基序列，人的很多模拟实验都要用到这些序列来研究特定排序的表达效果和功能，这时候不一定要从NCBI数据库下载真实的数据而只需要类似比例的序列进行模拟处理即可；

所以我就想到了这样一个小功能来帮助研究人员生成特定长度的特定数量的符合人类染色体一行碱基比例的DNA序列

#### getseqs函数

函数具体的实现思路很简单

```python
from random import *
pairwise = "AAATTTCCGG"
def getseq(num):
    seq = []
    for i in range(num):
        seq.append(pairwise[randint(0, 9)])
    return seq
```

因为查阅文献可以知道，一号染色体的AT与CG比例为3：2，所以我静态定义了一个字符串`pairwise`，其中的AT与CG比例正好符合要求，然后用生成随机数的方式从这个字符串里取出元素整合到新生成的序列中，然后执行这个函数规定次数，将所有生成的序列按行写入一个文件`seqs.txt`中

#### 使用手册

运行`smallprogram.py`文件，首先打开`On`，然后输入要求的长度和数量，点击上方按钮即可输出第一条符合要求的序列，然后在文件夹中的`seqs.txt`文件中即可找到其他所有满足该要求的序列

![](file:///home/egotist/PycharmProjects/dp_program/photo/4.jpg)

![](file:///home/egotist/PycharmProjects/dp_program/photo/5.jpg)

**如果想要生成特定比例的序列只需修改原文件中的`pairwise`字符串的AT与CG比例即可**

### 3.验证score服从的分布以及数据特征

#### 思路

这个是专业课的一个作业题，生成50条比例符合人类染色体1号序列的序列，然后任一两条之间进行比对生成score文件存储最终得分和seq文件存放最佳比对结果。

然后对于score文件绘制直方图判断其是否服从高斯分布，并计算score的平均值，标准差和进行k-s检验的结果



#### 具体实现

具体的实现函数和序列比对类似，只是多了`loadDta`函数读取`score.txt`文件并进行数据处理，`draw_hist`函数绘制直方图并存为`outcome.png`,然后将数据处理的结果全部输出到`outcome.txt`文件，并输出k-s检验的结果到软件中



#### 使用手册

运行`smallprogram.py`文件，首先打开`On`，然后点击中间上方按钮即可输出k-s检验结果，然后在文件夹中的`seq.txt,score.txt,outcome.txt,outcome.png`中可以查看详细结果

![](file:///home/egotist/PycharmProjects/dp_program/photo/6.jpg)

![](file:///home/egotist/PycharmProjects/dp_program/photo/9.jpg)

![](file:///home/egotist/PycharmProjects/dp_program/photo/10.jpg)

![](/home/egotist/PycharmProjects/dp_program/photo/7.jpg)





## 三.总结

实现的几个小功能其实都还比较简单，写起来没有啥问题.GUI编程对于所使用的包还不太熟悉，数据分析这一块还可以做很多东西，后续会逐渐补齐这一部分功能。

`tinker`包做这种可视化界面感觉终究是有一些麻烦的，并且美观性和功能性还有所欠缺，后续我会用`django`做一个网页的实现，然后加入更多使用的功能逐步完善这个小程序。

代码文件已上传至作者的github（egotist0），可以在里面获得实现文件，如果感兴趣或想提出一些建议的话,可以给本人的 Github 留言或直接参与修改
