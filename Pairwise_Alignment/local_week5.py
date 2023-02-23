#!/usr/bin/python3

# input files
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

#score function
def score(a,b):
	if a == "-" or b == "-":
		return -2
	elif  a != b:
		return -1
	elif a == b:
		return 1

# trace_matrix
def trace(score_matrix, trace_matrix):
	# if id == 0:
	# 	return 0
	# elif id == 1:
	# 	return 1
	# elif id == 2:
	# 	return 2
	# elif id == 3:
	# 	if ul == max(ul,l,u):
	# 		return 0
	# 	elif l == max(ul,l,u):
	# 		return 1
	# 	elif u == max(ul,l,u):
	# 		return 2
	i,j = len(dsq1)-1, len(dsq2)-1
	while i > 0 or j > 0:
		now = score_matrix[i][j]
		if i > 0 and j > 0:
			ul = score_matrix[i-1][j-1]
			l = score_matrix[i][j-1]
			u = score_matrix[i-1][j]
			score = max(ul,l,u)
			if ul == score:
				trace_matrix[i][j] = 0
				i = i-1
				j = j-1
			elif l == score:
				trace_matrix[i][j] = 1
				j = j-1
			elif u == score:
				trace_matrix[i][j] = 2
				i = i-1
		elif i > 0 and j <= 0:
			trace_matrix[i][j] = 2
			i = i-1
		elif i <= 0 and j > 0:
			trace_matrix[i][j] = 1
			j = j-1
	return trace_matrix


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

# print(' '.join(['%3s' % i for i in ' '+dsq2]))
# for i, p in enumerate(dsq1):
# 	line = [p] + [trace_matrix[i][j] for j in range(len(dsq2))]
# 	print(' '.join(['%3s' % i for i in line]))


#initialize matrix
#col for dsq1 and row for dsq2
score_matrix = {}
trace_matrix = {}

for i,p in enumerate(dsq1):
	score_matrix[i] = {}
	trace_matrix[i] = {}
	for j,q in enumerate(dsq2):
		if i == 0:
			score_matrix[i][j] = 0
			# trace: 1 for heng
			trace_matrix[i][j] = 1
			continue
		if j == 0:
			score_matrix[i][j] = 0
			# trace: 2 for shu
			trace_matrix[i][j] = 2
			continue
		ul = score_matrix[i-1][j-1] + score(p,q)
		# l for heng
		l = score_matrix[i][j-1] + score('-',q)
		# u for shu
		u = score_matrix[i-1][j]+score(p,'-')
		now = max([ul,l,u,0])
		score_matrix[i][j] = now
trace_matrix = trace(score_matrix,trace_matrix)
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
		middle = ' ' + middle
		dsq2 = dsq2[0: j-1]
	elif p == '2':
		align1 = dsq1[i-1] + align1
		align2 = '-' + align2
		middle = ' ' + middle
		dsq1 = dsq1[0: i-1]
similarity = numMatch / len(align1)

#print preety
i = 0
print('Local Alignment:')
print('seq 1: ' + str(title1))
print('seq 2: ' + str(title2))
print('Similarity: ' + str(similarity))
while i<len(align1):
	print (align1[i:i+97]+'\n'+middle[i:i+97]+'\n'+align2[i:i+97]+'\n')
	i = i+97

