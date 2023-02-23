#!/usr/bin/python3
import numpy as np
import time

T1 = time.clock()
# input files
# dsq2 must be reference sequence
fi2 = open("/home/linj2022/data/Hs_HBB_mRNA.fa", "r")
seq2 = fi2.readlines()
fi2.close()
title2 = seq2[0].strip("\n")
dsq2 = ""
i = 1
while (i < len(seq2)):
        dsq2 += seq2[i].strip("\n")
        i += 1
print(len(dsq2))

# dsq1 must bee query sequence
fi1 = open("/home/linj2022/data/Mm_HBB_mRNA.fa", "r")
seq1 = fi1.readlines()
fi1.close()
title1 = seq1[0].strip("\n")
dsq1 = "-"
i = 1
while (i < len(seq1)):
        dsq1 += seq1[i].strip("\n") 
        i += 1
print(len(dsq1))
dsq1_len = len(dsq1)-1

#score function
def score(a,b):
	if a == "-" or b == "-":
		return -2
	elif  a != b:
		return -1
	elif a == b:
		return 1

# trace_matrix
def trace(score_matrix):
	now = np.where(score_matrix == np.max(score_matrix))
	path_code = ""
	# for i in range(len(now)-1):
	i = len(now[0]) - 1
	x = now[0][i]
	y = now[1][i]
	bx, by = x,y
	while x > 0 and y > 0:
		ul = score_matrix[x-1][y-1]
		if ul == 0: break
		l = score_matrix[x][y-1]
		u = score_matrix[x-1][y]
		score = max(ul,l,u)
		if ul == score:
			x -= 1
			y -= 1
			path_code += "0"
		elif l == score:
			y -= 1
			path_code += "1"
		elif u == score:
			x -= 1
			path_code += "2"
	path_code = [bx, by, path_code]
	print(path_code)
	return path_code


# print(' '.join(['%3s' % i for i in ' '+dsq2]))
# for i, p in enumerate(dsq1):
# 	line = [p] + [trace_matrix[i][j] for j in range(len(dsq2))]
# 	print(' '.join(['%3s' % i for i in line]))


#initialize matrix
#col for dsq1 and row for dsq2
score_matrix = np.zeros([len(dsq1), len(dsq2)], int)
for i,p in enumerate(dsq1):
	for j,q in enumerate(dsq2):
		if i == 0 or j == 0: score_matrix[i][j] = 0
		else:	
			ul = score_matrix[i-1][j-1] + score(p,q)
			# l for heng
			l = score_matrix[i][j-1] + score('-',q)
			# u for shu
			u = score_matrix[i-1][j] + score(p,'-')
			now = max([ul,l,u,0])
			score_matrix[i][j] = now
path_code = trace(score_matrix)
# print matrix
# print(' '.join(['%3s' % i for i in ' '+dsq2]))
# for i, p in enumerate(dsq1):
# 	line = [p] + [score_matrix[i][j] for j in range(len(dsq2))]
# 	print(' '.join(['%3s' % i for i in line]))

#pathCode
#output:
#-ACGT
# | ||
#TA-GT
align1 = ''
middle = ''
numMatch = 0
align2 = ''
i = path_code[0]
j = path_code[1]
print(len(path_code[2]))
for p in path_code[2]:
	if p == '0':
		align1 = dsq1[i] + align1
		align2 = dsq2[j] + align2
		if dsq1[i-1] == dsq2[j-1]:
			middle = '|' + middle
			numMatch += 1
		else:
			middle = ' ' + middle
		dsq1 = dsq1[0: i]
		dsq2 = dsq2[0: j]
		i -= 1
		j -= 1
	elif p == '1':
		align1 = '-' + align1
		align2 = dsq2[j-1] + align2
		middle = ' ' + middle
		dsq2 = dsq2[0: j]
		j -= 1
	elif p == '2':
		align1 = dsq1[i-1] + align1
		align2 = '-' + align2
		middle = ' ' + middle
		dsq1 = dsq1[0: i]
		i -= 1
similarity = numMatch / dsq1_len

T2 = time.clock()
#print preety
i = 0
print('Local Alignment:')
print('seq 1: ' + str(title1))
print('seq 2: ' + str(title2))
print('Similarity: ' + str(similarity))
print("Time cost: " + str(T2-T1))
while i<len(align1):
	print (align1[i:i+97]+'\n'+middle[i:i+97]+'\n'+align2[i:i+97]+'\n')
	i = i+97

