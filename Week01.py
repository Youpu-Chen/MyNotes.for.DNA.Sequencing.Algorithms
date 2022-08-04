def readGenome(filename):
    genome = ''
    input = open(filename, 'r')
    for line in input:
        if not ">" in line:
            genome += line.strip()
    return genome

genome = readGenome('data/GCF_000840245.1_ViralProj14204_genomic.fna')
# print(genome)


def reverseComplement(seq):
    new_seq = ''
    dict1 = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for i in seq:
        new_seq += dict1[i]
    return new_seq[::-1]
# test
seq = 'AGGGGGGCCCCCCTATTTT'
# print(reverseComplement(seq))



def naive(p, t):
    '''this function is used to naive mathcing pattern in text'''
    occurrences = []
    
    for i in range(len(t) - len(p) + 1):
        match = True # Set original flag
        for j in range(len(p)):
            if not t[i + j] == p[j]:
                match = False
                break
        
        if match:
            occurrences.append(i)
    return occurrences

print('-------------------------------------------- Q1 --------------------------------------------')
'''Q1: How many times does AGGT or its reverse complement (ACCT) occur in the lambda virus genome?  
E.g. if AGGT occurs 10 times and ACCT occurs 12 times, you should report 22.'''
print(len(naive('AGGT', genome)) + len(naive('ACCT', genome)))


print('-------------------------------------------- Q2 --------------------------------------------')
'''Q2: How many times does TTAA or its reverse complement occur in the lambda virus genome? 
Hint: TTAA and its reverse complement are equal, so remember not to double count.'''
print(len(naive('TTAA', genome))) # TTAA and AATT have the same freq


print('-------------------------------------------- Q3 --------------------------------------------')
'''Q3: What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome?  
E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse complement ACTTAGT is at offset 29, then report 29.'''
# print(genome.find('ACTAAGT'))
# print(genome.find(reverseComplement('ACTAAGT')))
print('offset1: %d    offset2: %d' % (genome.find('ACTAAGT'), genome.find(reverseComplement('ACTAAGT'))))


print('-------------------------------------------- Q4 --------------------------------------------')
'''Q4: What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome?'''
# print(genome.find('AGTCGA'))
# print(genome.find(reverseComplement('AGTCGA')))
print('offset1: %d    offset2: %d' % (genome.find('AGTCGA'), genome.find(reverseComplement('AGTCGA'))))


print('-------------------------------------------- Q5 --------------------------------------------')
'''Q5: Rewrite the naive function to allow 2 base pair mismatches
naive-> naive_2mm'''
def naive_2mm(p, t):
    '''this function is used to naive mathcing pattern in text'''
    occurrences = []
    
    for i in range(len(t) - len(p) + 1):
        match = True # Set original flag
        count = 0
        for j in range(len(p)):
            if not t[i + j] == p[j]:
                count += 1
                if count > 2:
                    match = False
                
        if match:
            occurrences.append(i)
    return occurrences

print(len(naive('AGGT', genome)))
print(len(naive_2mm('AGGT', genome)))


print('-------------------------------------------- Q6 --------------------------------------------')
'''Q6: What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?'''
print(naive_2mm('AGGAGGTT', genome)[0])


print('-------------------------------------------- Q7 --------------------------------------------')
'''Q7: 
Finally, download and parse the provided FASTQ file containing real DNA sequencing reads derived from a human:

https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq

Note that the file has many reads in it and you should examine all of them together when answering this question.  The reads are taken from this study:

Ajay, S. S., Parker, S. C., Abaan, H. O., Fajardo, K. V. F., & Margulies, E. H. (2011). Accurate

and comprehensive sequencing of personal genomes. Genome research, 21(9), 1498-1505. 

This dataset has something wrong with it; one of the sequencing cycles is poor quality.

Report which sequencing cycle has the problem.  Remember that a sequencing cycle corresponds to a particular offset in all the reads. For example, if the leftmost read position seems to have a problem consistently across reads, report 0. If the fourth position from the left has the problem, report 3. Do whatever analysis you think is needed to identify the bad cycle. It might help to review the "Analyzing reads by position" video.'''

import matplotlib.pyplot as plt


def readfastq(filename):
    sequences = []
    qualities = []
    input = open(filename, 'r')

    while True:
        input.readline()
        seq = input.readline().strip()
        input.readline()
        qual = input.readline().strip()

        if len(seq) == 0:
            break
        sequences.append(seq)
        qualities.append(qual)
    input.close()
    return sequences, qualities

sequences, qualities = readfastq('data/ERR037900_1.first1000.fastq')
# print(sequences[:5])


def phred33ToQ(qual):
    '''Turn Phred+33 ASCII-encoded quality into Q'''
    return ord(qual) - 33
phred33ToQ('J')


def createHist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist
h = createHist(qualities)
plt.bar(range(len(h)), h)
plt.show()