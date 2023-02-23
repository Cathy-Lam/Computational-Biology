#!/usr/bin/python3

# input files
from re import A


fi1 = open("/home/linj2022/data/Hs_HBB_mRNA.fa", "r")
seq1 = fi1.readlines()
fi1.close()
title1 = seq1[0].strip("\n")
dsq1 = "-"
i = 1
while (i < len(seq1)):
        dsq1 += seq1[i].strip("\n")
        i += 1

fi2 = open("/home/linj2022/data/Mm_HBB_mRNA.fa", "r")
seq2 = fi2.readlines()
fi2.close()
title2 = seq2[0].strip("\n")
dsq2 = "-"
i = 1
while (i < len(seq2)):
        dsq2 += seq2[i].strip("\n")
        i += 1
# from itertools import count

# dsq1 = "-ATCGAGCTAGCGTA"
# dsq2 = "-TAGATTCTACGGTCAA"


#score function
def score(a,b):
	if a == "-" or b == "-":
		return -2 
	elif  a != b:
		return -1
	elif a == b:
		return 1


#traceback
def traceback(dsq1, dsq2, trace_matrix):
	i,j = len(dsq1)-1, len(dsq2)-1
	path_code = ''
	while i > 0 or j > 0:
		direction = trace_matrix[i][j]
		if direction == 0:
			i = i-1
			j = j-1
			# path_code = "0"+path_code
			path_code += "0"
		elif direction == 1:
			j = j-1
			# path_code = "1"+path_code
			path_code += "1"
		elif direction == 2:
			i = i-1
			# path_code = "2"+path_code
			path_code += "2"
	return path_code


#initialize matrix
#col for dsq1 and row for dsq2
score_matrix = {}
trace_matrix = {}

for i,p in enumerate(dsq1):
	# i for row and dsq1 in row
	score_matrix[i] = {}
	trace_matrix[i] = {}
	for j,q in enumerate(dsq2):
		# j for column and dsq2 in col
		if i == 0:
			score_matrix[i][j] = -j*2
			# trace: 1 for heng
			trace_matrix[i][j] = 1
			continue
		if j == 0:
			score_matrix[i][j] = -i*2
			# trace: 2 for shu
			trace_matrix[i][j] = 2
			continue
		ul = score_matrix[i-1][j-1] + score(p,q)
		# l for heng
		l = score_matrix[i][j-1] + score('-',q)
		# u for shu
		u = score_matrix[i-1][j]+score(p,'-')
		picked = max([ul,l,u])
		score_matrix[i][j] = picked
		trace_matrix[i][j] = [ul, l, u].index(picked)


# print matrix
# print(' '.join(['%3s' % i for i in ' '+dsq2]))
# for i, p in enumerate(dsq1):
# 	line = [p] + [score_matrix[i][j] for j in range(len(dsq2))]
# 	print(' '.join(['%3s' % i for i in line]))

# print(' '.join(['%3s' % i for i in ' '+dsq2]))
# for i, p in enumerate(dsq1):
# 	line = [p] + [trace_matrix[i][j] for j in range(len(dsq2))]
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
path_code = traceback(dsq1, dsq2, trace_matrix)
for p in path_code:
	i = len(dsq1)
	j = len(dsq2)
	if p == '0':
		align1 = dsq1[i-1] + align1
		align2 = dsq2[j-1] + align2
		if dsq1[i-1] == dsq2[j-1]:
			middle = '|' + middle
			numMatch += 1
		else:
			middle = ' ' + middle
		dsq1 = dsq1[0: i-1]
		dsq2 = dsq2[0: j-1]
	elif p == '1':
		align1 = '-' + align1
		align2 = dsq2[j-1] + align2
		middle = ' ' +  middle
		dsq2 = dsq2[0: j-1]
	elif p == '2':
		align1 = dsq1[i-1] + align1
		align2 = '-' + align2
		middle = ' ' + middle
		dsq1 = dsq1[0: i-1]
similarity = numMatch / len(align1)

#print preety
i = 0
print('Global Alignment:')
print('seq 1: ' + str(title1))
print('seq 2: ' + str(title2))
print('Similarity: ' + str(similarity))
while i<len(align1):
	print (align1[i:i+96]+'\n'+middle[i:i+96]+'\n'+align2[i:i+96]+'\n')
	i = i+96


