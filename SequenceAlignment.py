'''
Biological Sequence Alignment
@author: Gabriel A. Granell Jiménez
@date: November 4, 2022
Prof. Wilfredo E. Lugo
The Goal of project one is to implement your own Needleman-Wunsch algorithm using Dynamic Programming.
'''
from fileinput import filename
import sys
import csv
import numpy as np

#Using tabulation
def alignment(sequences,match_reward,mismatch_penalty,gap_penalty):
  if(len(sequences)>2):#If there are more than 2 sequences
    return

  if(len(sequences)==1):#If there is only one sequence.
    sequences.append("")

  if(len(sequences)==0):#If there are no sequences.
    print("0")
    return 

  #Initializing the two sequences stored.
  sequence_1=sequences[0]
  sequence_2=sequences[1]

  #Dimensions for matrix
  n=len(sequence_1)#Column dimension
  m=len(sequence_2) #Row dimension

  #Initializing the alignment matrices with zero.
  scoring_matrix=np.zeros((n+1,m+1))
  match_checker=np.zeros((n,m))

  #Filling the match checker matrix according to match or mismatch
  for i in range(n):
    for j in range(m):
        if sequence_1[i] == sequence_2[j]:
            match_checker[i][j]= match_reward #If they are the same, assign match reward.
        else:
            match_checker[i][j]= mismatch_penalty #If they are not the same, assign mismatch penalty.

  #Filling up the matrix using Needleman_Wunsch algorithm.
  for i in range(1,n+1):
    scoring_matrix[i][0]=i*gap_penalty

  for j in range(1,m+1):
    scoring_matrix[0][j]=j*gap_penalty

  #Modifying scoring matrix values
  for i in range(1,n+1):
    for j in range(1,m+1):
      scoring_matrix[i][j]=max(scoring_matrix[i-1][j-1] + match_checker[i-1][j-1],scoring_matrix[i][j - 1] + gap_penalty,scoring_matrix[i - 1][j] + gap_penalty)
      #First value in max is diagonal, second is horizontal, third is vertical.

  back_track(sequence_1,sequence_2,scoring_matrix,match_checker,n,m,gap_penalty)

#Back-tracking
def back_track(sequences1,sequences2,scoring_matrix,match_checker,n,m,gap_penalty):
  #Creating empty strings to store back-track
  aligned_1 = ""
  aligned_2 = ""

  #Assigning to new variable to not change value when printing
  x = n #Length of columns 
  y = m #Length of rows

  '''
  Special cases: When inputing through the calculator 
  in different ways these were some outputs the calculator 
  provided.They can be tested out by putting in csv:
    For 1) put: ",sequence(placeholder)"
    For 2) put: "sequence(placeholder),"
    For 3) put: ","
  '''
  if x==0 and y!=0:#1)If only an sequence 2 was inputed.
    for i in range(len(sequences2)):
      aligned_1+="-"
    
    aligned_2=sequences2
    print(aligned_1,aligned_2,int(scoring_matrix[x][y]))

  elif x!=0 and y==0:#2)If only sequence 1 was inputed.
    for i in range(len(sequences1)):
      aligned_2+="-"
    
    aligned_1=sequences1
    print(aligned_1,aligned_2,int(scoring_matrix[x][y]))

  elif x==0 and y==0:#3)If neither of them were inputed.
    print(aligned_1,aligned_2,int(scoring_matrix[x][y]))

  #Main back-tracking
  else:
    while(x >0 or y>0):
        #Directions
        diagonal=scoring_matrix[x-1][y-1]+ match_checker[x-1][y-1]
        left=scoring_matrix[x][y-1] + gap_penalty
        up=scoring_matrix[x-1][y] + gap_penalty

        if (diagonal>left and diagonal>up):#Will move diagonally
          aligned_1 = sequences1[x-1] + aligned_1
          aligned_2 = sequences2[y-1] + aligned_2
          x-=1
          y -= 1
        
        elif(left>=diagonal and left>=up):#Will move left (used >= to move like calculator given as example)
          aligned_1 = "-"+aligned_1
          aligned_2 = sequences2[y-1] + aligned_2

          y-=1

        else:
          aligned_1 =  sequences1[x-1] + aligned_1#Will move up
          aligned_2 = "-" + aligned_2

          x -= 1

    #Printing sequences and scoring of last row last column)
    print(aligned_1,aligned_2,int(scoring_matrix[n][m]))

def main():
    if len(sys.argv)>1:
      filename=sys.argv[1]
      file=open(sys.argv[1])
      reader=csv.reader(file)
      next(reader,None)
      for i in reader:
        alignment(i,1,-1,-2)#sequences,match_reward,mismatch_penalty,gap_penalty

if __name__ == '__main__':
  main()