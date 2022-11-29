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

- Biobase sequences consist of ATCG, a scoring rule needs to be set for comparison, I set `match=5,dismatch=-4,gap=-2.5` here, if the user needs to modify the rule, the location of the modified rule is in lines 26-28 of `data_analysis.py` file and lines 16-18 of `dp.py` file. Lines 16-18 of the `dp.py` file. ![](https://github.com/egotist0/Gene_factory/tree/master/photo/1.jpg)

![](https://github.com/egotist0/Gene_factory/tree/master/photo/2.jpg)



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

#### backtracking function

+ the way backtracking works is to look at the top left, top and left maximum positions of each backtracking locus, and finally the entire backtracking path is obtained.

+ The `printAlign` function has more parameters, `score` is the input two-dimensional scoring array, `i` and `j` represent the position of the scoring array backtracked at this time, the initial position is from (len1,len2), `s1` and `s2` are the two sequences of backtracking, `saln` and `raln` store the final output Best match result

+ The three `if` judgments are to determine which score is highest on the top left, left, or top, and the highest is the optimal matching path: 1.

  1. starting from the bottom right, if the maximum value appears on top, a GAP ("-") is introduced for this sequence horizontally, and the bases are taken there for this sequence vertically.
  2. if the maximum value appears on the left, introduce a GAP ("-") into this sequence vertically and take the base there for this sequence horizontally; 
  3. if the maximum value appears in the upper left corner, no GAP is introduced, and the bases are taken both vertically and horizontally.
  4. Of course, the path of this backtracking is not necessarily unique, when two of the three positions are the same, both paths are feasible, but for the convenience of the output of the applet, I use `if..... .elif` so that only one sequence will be sifted out

#### manual

Run `smallprogram.py` file, first open `On`, then input two sequences, click the button above to output the highest score and the best match

![](https://github.com/egotist0/Gene_factory/tree/master/photo/photo/%E5%88%9D%E5%A7%8B.jpg)

![](https://github.com/egotist0/Gene_factory/tree/master/photo/3.png)

**Since this is an excellent text string matching algorithm, it can be used not only for base sequence matching, but similar text matching, similarity matching, text checking, protein sequence matching can be performed as well**



### 2. Generate random sequences scaled to human chromosomes

#### Idea

(a) Human chromosome 1 has 230 million base sequences, which are used in many human simulation experiments to study the expression effects and functions of specific sequences, when it is not necessary to download real data from the NCBI database but only a similar proportion of sequences for simulation.

So I came up with this little function to help researchers generate a specific number of DNA sequences of a specific length that match the ratio of bases in a row of human chromosomes

#### getseqs function

The exact implementation of the function is simple

```python
from random import *
pairwise = "AAATTTCCGG"
def getseq(num):
    seq = []
    for i in range(num):
        seq.append(pairwise[randint(0, 9)])
    return seq
```

Since a review of the literature shows that the AT to CG ratio of chromosome 1 is 3:2, I statically define a string `pairwise` in which the AT to CG ratio is exactly as required, and then use the generation of random numbers to integrate the elements from this string into the newly generated sequence, and then execute this function a prescribed number of times to write all the generated sequences by line to a file `seqs.txt`.

#### User Manual

Run the `smallprogram.py` file, first open `On`, then enter the required length and number, click the button above to output the first sequence that meets the requirement, then you can find all other sequences that meet the requirement in the `seqs.txt` file in the folder

![](https://github.com/egotist0/Gene_factory/tree/master/photo/4.jpg)

![](https://github.com/egotist0/Gene_factory/tree/master/photo/5.jpg)

**If you want to generate a sequence with a specific ratio just modify the AT to CG ratio of the `pairwise` string in the original file**

### 3. Verify the distribution of score obedience and the data characteristics

#### Ideas

This is a homework problem in a professional course. Generate 50 sequences that match the human chromosome 1 sequence in proportion, and then compare between any two of them to generate a score file to store the final score and a seq file to store the best comparison results.

Then the histogram is plotted for the score file to determine whether it follows a Gaussian distribution, and the mean, standard deviation and k-s test results of the score are calculated.



#### Concrete implementation

The specific implementation functions are similar to the sequence matching, except that the `loadDta` function reads the `score.txt` file and processes the data, the `draw_hist` function draws the histogram and saves it as `outcome.png`, then outputs all the results of data processing to the `outcome.txt` file, and outputs the results of the k-s test to the software



#### User Manual

Run the `smallprogram.py` file, first open `On`, then click the top middle button to output the k-s test results, then you can view the detailed results in the folder `seq.txt,score.txt,outcome.txt,outcome.png`

![](https://github.com/egotist0/Gene_factory/tree/master/photo/6.jpg)

![](https://github.com/egotist0/Gene_factory/tree/master/photo/9.jpg)

![](https://github.com/egotist0/Gene_factory/tree/master/photo/10.jpg)

![](https://github.com/egotist0/Gene_factory/tree/master/photo/7.jpg)





## III. Summary

The GUI programming is not too familiar with the package used, and there is still a lot that can be done for data analysis, so I will gradually make up for this part of the functionality.

I will use `tinker` package to do this kind of visualization interface feel that there is some trouble after all, and the beauty and functionality is still lacking, I will use `django` to do a web page implementation, and then add more features to gradually improve this small program.

The code file has been uploaded to the author's github (egotist0), you can get the implementation file inside, if you are interested or want to make some suggestions, you can leave a message to my Github or directly participate in the modification
