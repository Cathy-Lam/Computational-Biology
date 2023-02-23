#!/usr/bin/python3
import time
import numpy as np

# it's the version with hash HBB and genome

# MaxPQ class
# initialize: maxPQ = MaxPQ()
# then use insert: maxPQ.insert(int)
# or insert list: maxPQ.insertList(list)
# only mantain the top TEN queue
class MaxPQ:
    def __init__(self):
        self.size = 0
        self.pq = [0]

    def max(self):
        if len(self.pq) > 1:
            return self.pq[1]
        else: return -1
    
    def swap(self, i:int, j: int):
        temp = self.pq[i]
        self.pq[i], self.pq[j] = self.pq[j], temp
    
    def less(self, i:int, j:int) ->bool:
        return self.pq[i] < self.pq[j]
    
    def left(self, i:int) ->int:
        return i * 2

    def right(self, i:int) ->int:
        return i * 2 + 1

    def parent(self, i:int) ->int:
        return i // 2
    
    def swim(self, x:int):
        f = self.parent(x)
        while x > 1 and self.less(f, x):
            self.swap(f, x)
            x = f
            f = self.parent(x)
    
    def sink(self, x:int):
        while self.left(x) <= self.size:
            m = self.left(x)
            r = self.right(x)

            if r <= self.size and self.less(m, r):
                m = r
            if self.less(m, x):
                break
            
            self.swap(x, m)
            x = m
    
    def insert(self, x):
        self.size += 1
        self.pq.append(x)
        self.swim(self.size)
    
    def delMax(self):
        m = self.pq[1]
        self.swap(1, self.size)
        del self.pq[self.size]
        self.size -= 1
        self.sink(1)
        return m
    
    def delMin(self):
        del self.pq[self.size]
        self.size -= 1

    def insertList(self, insertList:list):
        n = len(insertList)
        i = 0
        while i < n:
            self.insert(insertList[i])
            i += 1

    def topNumInsert(self, num):
        while self.size > num: self.delMin()

    def topNumInsertList(self, insertList:list, num):
        i = num + 1
        if len(insertList) <= num:
            self.insertList(insertList)
        else:
            self.insertList(insertList)
            self.topNumInsert(num)

    def print(self):
        return sorted([i for i in self.pq[1:]], reverse=True)

#binary search
def binary_search_left(value, list):
	lo = 0
	hi = len(list) - 1
	while lo < hi:
		mid = int(lo + (hi - lo + 1)/2)
		if list[mid] <= value:
			lo = mid
		else:
			hi = mid - 1
	return lo if list[lo] <= value else -1

def binary_search_right(value, list):
	lo = 0
	hi = len(list) - 1
	while lo < hi:
		mid = int(lo + (hi - lo)/2)
		if list[mid] >= value:
			hi = mid
		else:
			lo = mid + 1
	return lo if list[lo] >= value else -1

# do hash with HBB and HBB_inverse
# key: seq, value: index
# the hash is too big, it must need some approach to improve it
# do we really need a only one-bp-moving splits
def doHash(seq):
	seq_hash = {}
	i = 0
	while (i+11) <= len(seq):
		# if maskLowComplexity:
		if seq[i:i+11] in seq_hash:
			seq_hash[seq[i:i+11]].append(i)
		else:
			seq_hash[seq[i:i+11]] = [i]
		i += 1
	seq_hash[seq[i:len(seq)]] = [i]
	return seq_hash

# do compare to get primary seed
# principle: same
# it need to score the splits
# format: seed[chr] = [sorte whole match genome_id]
# maybe arrange as MaxPQ? or just sort like shell or binary since we just need the top ten seq
def doSeed(seq, genome_hash):
	seed = []
	k = 0
	while (k + 11) <= len(seq):
		sq = seq[k:k+11]
		if sq in genome_hash:
			seed += genome_hash[sq] # genome_hash[sq] = [id]
		k += 1
	seed.sort()
	return seed

# filter
def findNum(seed:list, index, up:int) -> list:
	m = 1
	j = 0
	while index + j < len(seed):
		if seed[index + j] <= up: j += 1
		else: 
			break
	m += j
	return [m, seed[index], seed[min((index+j),len(seed)-1)]]

# format: seed_filtered[m] = [[0/1, chr, id], [0/1, chr, id]]
# seed[chr] = [match list]
# 0 for HBB, 1 for HBB_inverse
# build a MaxPQ tree to maintain the first 100 matchNum id
def doSeedFilter(seq_len, seed, seed_inverse, genome):
	# seed_primary = {}
	# seed_filtered = {}
	# match_MaxPQ = MaxPQ()
	# for chr in seed: #0
	# 	j = 0
	# 	matchList = [0, 0, 0] #[m, genome_id]
	# 	temp_j = matchList[1]
	# 	# while j < len(seed[chr]):
	# 	for j,id in enumerate(seed[chr]):
	# 		up = min(id + (1.5 * seq_len), len(genome[chr]))
	# 		matchList = findNum(seed[chr], j, up)
	# 		m = matchList[0]
	# 		if m > 32 and (temp_j+3) < matchList[2]:
	# 			temp_j = matchList[2]
	# 			if m in seed_primary: seed_primary[m].append([0, chr, matchList[1]])
	# 			else: seed_primary[m] = [[0, chr, matchList[1]]]
	# 		if up == len(genome[chr]): break
	

	# for chr in seed_inverse: #1
	# 	matchList = [0, 0, 0] #[m, id]	
	# 	temp_j = matchList[1]
	# 	j = 0
	# 	for j,id in enumerate(seed_inverse[chr]):
	# 		up = min(id + (1.5 * seq_len), len(genome[chr]))
	# 		matchList = findNum(seed_inverse[chr], j, up)
	# 		m = matchList[0]
	# 		if m > 32 and (temp_j+3) < matchList[2]:
	# 			temp_j = matchList[2]
	# 			if m in seed_primary: seed_primary[m].append([1, chr, matchList[1]])
	# 			else: seed_primary[m] = [[1, chr, matchList[1]]]
	# 		if up == len(genome[chr]): break

	# numInsert = 0
	# for i, m in enumerate(seed_primary):
	# 	numInsert += len(seed_primary[m])
	# 	if numInsert < 100: 
	# 		match_MaxPQ.insert(m)
	# 	else: match_MaxPQ.topNumInsert(m, numInsert)
	
	# matchList = {}
	# matchList_inverse = {}
	# # matchList_RD = {}
	# seed_primary = {}
	# seed_filtered = {}
	# match_MaxPQ = MaxPQ()
	# for chr in seed: #0
	# 	matchList[chr] = []
	# 	# matchList_RD[chr] = []
	# 	# matchList = [0, 0, 0] #[m, genome_id, genome_end_id]
	# 	for j,id in enumerate(seed[chr]):
	# 		up = min(id + (1.5 * seq_len), len(genome[chr]))
	# 		matchList[chr].append(findNum(seed[chr], j, up)) #[m, genome_id, genome_end_id]
	# 		# matchList_RD[chr].append(matchList[0])

	# for chr in seed_inverse: #0
	# 	matchList_inverse[chr] = []
	# 	# matchList_RD[chr] = []
	# 	# matchList = [0, 0, 0] #[m, genome_id, genome_end_id]
	# 	for j,id in enumerate(seed_inverse[chr]):
	# 		up = min(id + (1.5 * seq_len), len(genome[chr]))
	# 		matchList_inverse[chr].append(findNum(seed_inverse[chr], j, up)) #[m, genome_id, genome_end_id]
	# 		# matchList_RD[chr].append(matchList[0])
	# 	# 	if m > 32 and (temp_j+3) < matchList[2]:
	# 	# 		if m in seed_primary: seed_primary[m].append([0, chr, matchList[1]])
	# 	# 		else: seed_primary[m] = [[0, chr, matchList[1]]]
	# 	# 	if up == len(genome[chr]): break
	
	# for chr in matchList: # take the largest
	# 	for i,id in enumerate(matchList[chr]):
	# 		# left = max(i,0)
	# 		# right = min(i,len(matchList_inverse[chr])-1)
	# 		if i == 0 and i == len(matchList[chr])-1: continue 
	# 		if id[0] > matchList[chr][i-1][0] and id[0] >matchList[chr][i+1][0] and id[0]>32:
	# 			if id[0] in seed_primary: seed_primary[id[0]].append([0, chr, id[1]])
	# 			else: seed_primary[id[0]] = [[0, chr, id[1]]]

	# for chr in matchList_inverse: # take the largest
	# 	for i,id in enumerate(matchList_inverse[chr]):
	# 		# left = max(i,0)
	# 		# right = min(i,len(matchList_inverse[chr])-1)
	# 		if i == 0 and i == len(matchList_inverse[chr])-1: continue
	# 		if id[0] > matchList_inverse[chr][i-1][0] and id[0] > matchList_inverse[chr][i+1][0] and id[0]>32:
	# 			if id[0] in seed_primary: seed_primary[id[0]].append([1, chr, id[1]])
	# 			else: seed_primary[id[0]] = [[1, chr, id[1]]]


	seed_primary = {}
	seed_filtered = {}
	match_MaxPQ = MaxPQ()
	for chr in seed: #0
		temp_List_left = []
		temp_List_middle = []
		temp_List_right = []
		for j,id in enumerate(seed[chr]):
			up = min(id + (1.5 * seq_len), len(genome[chr]))
			temp_List_left = temp_List_middle
			temp_List_middle = temp_List_right
			temp_List_right = findNum(seed[chr], j, up) #[m, genome_id, genome_end_id]
			if temp_List_left == []: continue
			elif temp_List_middle[0] > temp_List_left[0] and temp_List_middle[0] > temp_List_right[0] and temp_List_middle[0] > 32:
				if temp_List_middle[0] in seed_primary: seed_primary[temp_List_middle[0]].append([0, chr, temp_List_middle[1]])
				else: seed_primary[temp_List_middle[0]] = [[0, chr, temp_List_middle[1]]]


	for chr in seed_inverse: #0
		temp_List_left = []
		temp_List_middle = []
		temp_List_right = []
		# matchList_RD[chr] = []
		# matchList = [0, 0, 0] #[m, genome_id, genome_end_id]
		for j,id in enumerate(seed_inverse[chr]):
			up = min(id + (1.5 * seq_len), len(genome[chr]))
			temp_List_left = temp_List_middle
			temp_List_middle = temp_List_right
			temp_List_right = findNum(seed_inverse[chr], j, up) #[m, genome_id, genome_end_id]
			if temp_List_left == []: continue
			elif temp_List_middle[0] > temp_List_left[0] and temp_List_middle[0] > temp_List_right[0] and temp_List_middle[0] > 32:
				if temp_List_middle[0] in seed_primary: seed_primary[temp_List_middle[0]].append([1, chr, temp_List_middle[1]])
				else: seed_primary[temp_List_middle[0]] = [[1, chr, temp_List_middle[1]]]

	numInsert = 0
	match_MaxPQ.insertList([m for i, m in enumerate(seed_primary)])
	MaxList = match_MaxPQ.print()
	for i, m in enumerate(MaxList):
		if m in seed_primary:
			numInsert += len(seed_primary[m])
			if numInsert >= 30: break # filter index
			seed_filtered[m] = seed_primary[m]
	return seed_filtered

# brute force: local alignment: all match>4 seq

#score function
def localAlignScore(a,b):
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
	return path_code

# single local align
# can change the score method: (now)only similarity
# format: alignResult = [position, align1_correct, middle_correct, align2_correct, similarity]
# initialize at last
def localAlign(id, dsq1, dsq2):
	dsq2_len = len(dsq2) - 1
	score_matrix = np.zeros([len(dsq1), len(dsq2)], int)
	for i,p in enumerate(dsq1):
		for j,q in enumerate(dsq2):
			if i == 0 or j == 0: score_matrix[i][j] = 0
			else:	
				ul = score_matrix[i-1][j-1] + localAlignScore(p,q)
				# l for heng
				l = score_matrix[i][j-1] + localAlignScore('-',q)
				# u for shu
				u = score_matrix[i-1][j] + localAlignScore(p,'-')
				now = max([ul,l,u,0])
				score_matrix[i][j] = now
	path_code = trace(score_matrix)

	align1 = ''
	middle = ''
	align2 = ''
	i = path_code[0]
	j = path_code[1]
	moveInGenome = 0
	numMatch = 0
	for p in path_code[2]:
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
			moveInGenome += 1
			i -= 1
			j -= 1
		elif p == '1':
			align1 = '-' + align1
			align2 = dsq2[j-1] + align2
			middle = ' ' + middle
			dsq2 = dsq2[0: j-1]
			moveInGenome += 1
			j -= 1
		elif p == '2':
			align1 = dsq1[i-1] + align1
			align2 = '-' + align2
			middle = ' ' + middle
			dsq1 = dsq1[0: i-1]
			i -= 1
	position_begin = id + path_code[1] - moveInGenome
	position_end = position_begin + moveInGenome
	position = str(position_begin) + " - " + str(position_end)
	similarity = numMatch / dsq2_len
	alignResult = [position, align1, middle, align2, similarity]
	return alignResult

# whole local align sum with input seq_match
# final score: E = Kmn e^(-λS)
# E = the expect value = the number of highscoring segment pairs (HSPs) expected
# m, n = the length of two sequences
# λ, K = Karlin Altschul statistics
# initialize at last
# try set(): a & b to extend the seed and use the length to sort
# now: use the number of match to sort
# format: seed_filtered[m] = [[0/1, chr, id], [0/1, chr, id]]
# format: alignResult[similarity] = [chr, index(id+position), align1_correct, middle_correct, align2_correct]
def doLocalAlign(seq, seq_inverse, seed, genome):
	seq_len = len(seq)
	alignResult = {}
	for i,m in enumerate(seed):
		j = 0
		while j < len(seed[m]):
			mList = seed[m][j]
			if mList[0] == 0: dsq2 = '-' + seq
			else: dsq2 = '-' + seq_inverse
			up = min(mList[2] + (1.5 * seq_len), len(genome[mList[1]]))
			dsq1 = '-' + genome[mList[1]][mList[2]:int(up)]
			tempAlignResult = localAlign(mList[2], dsq1, dsq2)
			if tempAlignResult[-1] in alignResult:
				alignResult[tempAlignResult[-1]].append([mList[1]] + tempAlignResult[0:4])
			else: 
				alignResult[tempAlignResult[-1]] = [[mList[1]] + tempAlignResult[0:4]]
			j += 1
	
	# initialize
	dsq1, dsq2 = '', ''

	return alignResult

# print format
def printResult(MaxList, MaxSimilarity):
	print (str(MaxList[0]) + " " + str(MaxList[1]) + "   " + "similarity: " + str(MaxSimilarity))
	j = 0
	while j < len(MaxList[2]):
		print (MaxList[2][j:j+96] + '\n' + MaxList[3][j:j+96] + '\n' + MaxList[4][j:j+96]+'\n')
		j += 96

# use max priority queue to sort the number of match
# format: alignResult[similarity][chr, id, align1, middle, align3]
# maxPQ only about the similarity
def maxPQ_score(alignResult):
	similarity_MaxPQ = MaxPQ()
	insertLiist = [similarity for i,similarity in enumerate(alignResult)]
	similarity_MaxPQ.topNumInsertList(insertLiist, 10)

	similarityTopList = similarity_MaxPQ.print()
	for i, MaxSimilarity in enumerate(similarityTopList):
		MaxSimilarity = similarityTopList[i]
		if MaxSimilarity in alignResult:
			if len(alignResult[MaxSimilarity]) == 1:
				printResult(alignResult[MaxSimilarity][0], MaxSimilarity)
			else:
				for lst in alignResult[MaxSimilarity]:
					printResult(lst, MaxSimilarity)

T1 = time.clock()

# input HBB sequence
fi = open("/home/linj2022/data/HBB_p13.fa", "r")
HBB_temp = fi.readlines()
fi.close()
HBB_title = HBB_temp[0].strip("\n")
HBB = ""
hbb_i = 1
while (hbb_i < len(HBB_temp)):
	HBB += HBB_temp[hbb_i].strip("\n")
	hbb_i += 1

# get the inverse HBB
complementaryBase = {'A':'T','T':'A','G':'C','C':'G'}
HBB_inverse = ""
for base in HBB:
	HBB_inverse += complementaryBase[base]
HBB_inverse = HBB_inverse[::-1]

end = False
genome = {}
genome_hash = {}
title = ""
seq = ""
HBB_len = len(HBB)
HBB_seed = {}
HBB_inverse_seed = {}
with open("/data/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa", "r") as fi:
	for line in fi:
		if ">chr" in line:
			if title == "":
				title = line.strip("\n")

			if seq != "":
				# is it the end chr
				if "Y" in title:
					end = True

				genome[title] = seq
				genome_hash[title] = doHash(genome[title])
				HBB_seed[title] = doSeed(HBB, genome_hash[title])
				HBB_inverse_seed[title] = doSeed(HBB_inverse, genome_hash[title])
				genome_hash[title] = []

				# initialize seq
				seq = ""
				title = line.strip("\n")
				
		else:
			if end == False:
				seq += line.strip("\n")
			else:
				if (">" in line):
					break
				else:
					seq += line.strip("\n")
genome[title] = seq
genome_hash[title] = doHash(genome[title])
HBB_seed[title] = doSeed(HBB, genome_hash[title])
HBB_inverse_seed[title] = doSeed(HBB_inverse, genome_hash[title])
genome_hash[title] = []
seq = ""

T2 = time.clock()
print("time cost for read and doHash and doSeed: " + str(T2-T1))

seed_filtered = doSeedFilter(HBB_len, HBB_seed, HBB_inverse_seed, genome)

alignResult = doLocalAlign(HBB, HBB_inverse, seed_filtered, genome) #seq, seq_inverse, seed, genome
seed_filtered = {}
HBB_inverse_seed = {}
genome = {}
T3 = time.clock()
print("time cost for doLocalAlign with filter: " + str(T3-T2))

maxPQ_score(alignResult)
T4 = time.clock()
print("time cost for final maxPQ and print top TEN: " + str(T4-T3))
print("whole time cost: " + str(T4 - T1))