#!/usr/bin/python3
end = False
end21 = True
genome = {}
title = ""
seq = ""
with open("/data/software/refdata-gex-GRCh38-2020-A/fasta/genome.fa", "r") as fi:
	for line in fi:
		if "chr11" in line:
			end = True
			title = line.strip("\n")
				
		else:
			if end == False:
				continue
			else:
				if (">" in line):
					break
				else:
					seq += line.strip("\n")
ge = seq.strip("N")
print(str(title))
print(str(ge))

