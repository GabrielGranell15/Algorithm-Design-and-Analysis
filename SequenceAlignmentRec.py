'''
Biological Sequence Alignment
@author: Gabriel A. Granell Jiménez
@date: November 4, 2022
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
  

  aligned_1 = ""
  aligned_2 = ""
  #Special cases: For example if only 1 sequence was inputed or none
  if n==0 and m!=0:
    for i in range(len(sequence_2)):
      aligned_1+="-"
    aligned_2=sequence_2
    print(aligned_1,aligned_2,int(scoring_matrix[n][m]))

  elif n!=0 and m==0:#2)If only sequence 1 was inputed.
    for i in range(len(sequence_1)):
      aligned_2+="-"
    
    aligned_1=sequence_1
    print(aligned_1,aligned_2,int(scoring_matrix[n][m]))

  elif n==0 and m==0:#3)If neither of them were inputed.
    print(aligned_1,aligned_2,int(scoring_matrix[n][m]))
  
  else:#Main backtracking
    arr=backtracking(aligned_1,aligned_2,sequence_1,sequence_2,scoring_matrix,match_checker,n,m,gap_penalty)
    print(arr[0],arr[1],int(scoring_matrix[n][m]))

#Recursive backtracking
def backtracking(aligned_1,aligned_2,sequence1,sequence2,scoring_matrix,match_checker,n,m,gap_penalty):
  
  if (0<n or 0<m):
    diagonal=scoring_matrix[n-1][m-1]+ match_checker[n-1][m-1]
    left=scoring_matrix[n][m-1] + gap_penalty
    up=scoring_matrix[n-1][m] + gap_penalty

    if(diagonal>left and diagonal>up):#Will move diagonally
      return backtracking(sequence1[n-1] + aligned_1,sequence2[m-1] + aligned_2,sequence1,sequence2,scoring_matrix,match_checker,n-1,m-1,gap_penalty)

    elif(left>=diagonal and left>=up):#Will move left (used >= to move like calculator given as example)
      return backtracking("-"+aligned_1,sequence2[m-1] + aligned_2,sequence1,sequence2,scoring_matrix,match_checker,n,m-1,gap_penalty)
    
    else:#Will move up
      return backtracking(sequence1[n-1] + aligned_1,"-" + aligned_2,sequence1,sequence2,scoring_matrix,match_checker,n-1,m,gap_penalty)
  
  else:#Returning arr with aligned sequences
    return aligned_1,aligned_2
    

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