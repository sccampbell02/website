---
date: "2019-12-14"
math: "true"
title: Sequence Alignment with Python
---

Nucleotide sequence alignment can be tedious to code by hand, but thanks to services like [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&BLAST_SPEC=blast2seq&LINK_LOC=align2seq) and the Python package biopython, understanding the alignment process is not necessary for producing parsimonious alignments. Even though hand-coding isn't the quickest method, understanding traditional algorithms can clear up how modern programs work and aid in troubleshooting. I'll go over a common global alignment algorithm, the Needleman-Wunsch algorithm.

## Needleman-Wunsch algorithm
This is a global alignment algorithm, identifying the best alignment after considering each residue within the entire sequence. It is done by completing an alignment matrix of two sequences according to an arbitrary scoring system. I'll walk through an example of doing this on paper.

Consider the sequences **ATTACCA** and **ATACGA**. Add a "gap" character (-) to the front of each sequence, and assign each to the first row or column of an alignment matrix:

{{% portfolio image="/images/Presentation1.png"%}}

{{% /portfolio %}}

We will esentially be assigning a numerical score to each of the remaining empty cells that depends on whether the corresponding nucleotides match up. The score will be calculated based on whether there is a match, mismatch, or gap (leaving a gap in one sequence where there is a nucleotide in the other; this will make more sense later on). Let's use a scoring system where a match = +1, mismatch = -1, and gap penalty = -1.

Start by filling in the intersection of the two gaps with a 0. Then, complete the second row from left to right. Because we are moving horizontally (i.e., matching all values to a gap), these cells all incur the gap penalty. Thus, the value of each cell is 1 less than the cell to its immediate left. Draw an arrow for each cell pointing to the cell from which its value was calculated.

{{% portfolio image="/images/SW2.gif"%}}

{{% /portfolio %}}

Because the gap penalty is also incurred for vertical movement, repeat this process for the second column moving from top to bottom.

{{% portfolio image="/images/Slide7.png"%}}

{{% /portfolio %}}

There are three possible values for each of the remainder of empty cells. This value could be the value of the cell above minus one (vertical movement, so gap penalty), value of the upper left cell plus or minus one, depending on whether the focal cell is an intersection of matching or mismatching nucleotides, or the value of the cell to the left minus one (horizontal movement, so gap penalty). The highest of these three values is recorded.

Consider the cell in row 3, column 3. This is the intersection of two A's: a match. The cell above is -1, so -2 (-1-(2)) is a possible value. The cell above and to the left is 0, and because our focal cell is a match 1 is a possible value. The cell to the left is -1, so -2 is once again a possible value. 1 is greater than -2, so 1 is the focal cell's value. Draw an arrow from the focal cell to the cell that produced our final value (top left).

{{% portfolio image="/images/SW3.gif"%}}

{{% /portfolio %}}

The process is repeated for the cell in the third column, fourth row. This is the intersection of an A and a T: a mismatch. The cell above is 1, so 0 (1-1) is a possible value. The cell above and to the left is -2, so -3 (-2-(1)) is a possible value due to the mismatch. The cell to the left is -3, so -4 is a possible value. 0 is the greatest of all these values, so we record 0 and draw an arrow from the focal cell to the one above.

{{% portfolio image="/images/SW4.gif"%}}

{{% /portfolio %}}

This process is repeated until the matrix is complete. If two surrounding cells both produce the same focal cell value, and this value is the highest of all possible values, two arrows are drawn.

{{% portfolio image="/images/SW5.gif"%}}

{{% /portfolio %}}

Once all arrows are drawn, start at the bottom right cell and follow the arrows, tracing the path. In this case, the path branches at one point.

{{% portfolio image="/images/SW6.gif"%}}

{{% /portfolio %}}

Skipping the gap intersection, move from top left to bottom right. Whatever two nucleotides intersect along the path are aligned. If the path goes along the horizontal for more than one cell, a gap is inserted in the vertical sequence, and vice versa for vertical movement. There are two possible paths and thus two equally viable alignments.

{{% portfolio image="/images/SW7.gif"%}}

{{% /portfolio %}}

Now to translate this to a Python function. This algorithm requires us to determine whether an intersection is a match, which could be aided by  a match helper function:

``` python
def match_fun(a,b): # Helper function to define match/mismatch scores
 if a==b:
 score=1 #match score
 else:
 score=-1 #mismatch score
 return score
```
Now, to create the main function. As on paper, we'l start by creating an empty array with a length one unit longer than one sequence, and a width one unit longer than the other sequence.

``` python
def needleman(seq1,seq2,gap_penalty=-1): # Needleman-Wunsch algorithm in a funct
  m=len(seq1) #length of horizontal sequence
  n=len(seq2) #length of vertical sequence
  score=np.empty(shape=[n+1,m+1]) #array to hold scores
```

We'll fill in the initial gap penalty squares, then iterate through the remainder of the cells, calculating the three possible values and settling on the maximum.

``` python
def needleman(seq1,seq2,gap_penalty=-1): # Needleman-Wunsch algorithm in a function
  m=len(seq1) #length of horizontal sequence
  n=len(seq2) #length of vertical sequence
  score=np.empty(shape=[n+1,m+1]) #array to hold scores
  for j in range(0, m + 1): score[0][j] = gap_penalty * j
  for i in range(0, n + 1): score[i][0] = gap_penalty * i
  for i in range(1, n + 1):
    for j in range(1, m + 1):
      insert = score[i - 1][j] + gap_penalty
      delete = score[i][j - 1] + gap_penalty
      match = score[i - 1][j - 1] + match_fun(seq1[j-1], seq2[i-1])
      score[i][j] = max(match, delete, insert)
  return score
```
This returns a matrix much like the one we filled out by hand. For the alignment, we'll create one final function that takes this matrix as an input. Starting at the bottom right hand corner, an if/else statement determines whether the value was determined by a match/mismatch to the diagonal or a gap penalty to the top or left. We then iterate to the chosen cell and repeat. Diagonal movement records an alignment between the sequences, while vertical or horizontal movement inserts a gap in one of the two sequences.
``` python
def get_alignment(seq1,seq2,mat,gap_penalty=-1):
    j=len(seq1) #length of horizontal sequence
    i=len(seq2) #length of vertical sequence
    AlignA=""   #empty strings to hold alignments
    AlignB=""
    while i > 0 or j > 0:  #start at the bottom right corner: if from diagonal, align sequences
        if i > 0 and j > 0 and mat[i,j]==mat[i-1][j-1]+match_fun(seq1[j-1],seq2[i-1]):
            AlignA = seq1[j-1] + AlignA
            AlignB = seq2[i-1] + AlignB
            i -= 1
            j -= 1
        elif j > 0 and mat[i,j]==mat[i][j-1]+gap_penalty: #if from above, put gap in vertical sequence
            AlignA = seq1[j-1] + AlignA
            AlignB = "-" + AlignB
            j -= 1
        else: #if from the left, put gap in horizontal sequence
            AlignA = "-" + AlignA
            AlignB = seq2[i-1]+AlignB
            i -= 1
    return AlignA, AlignB #return both alignments
```
Let's try this on two example sequences:

``` python
s1="ACGCTTACCG"
s2="AGCCTACCCC"
get_alignment(s1,s2,needleman(s1,s2))
```
**## ('ACGCTTA-CCG', 'A-GCCTACCCC')**
