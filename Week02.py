"""bm_preproc.py: Boyer-Moore preprocessing."""

__author__ = "Ben Langmead"

import unittest


def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]


class TestBoyerMoorePreproc(unittest.TestCase):

    def test_z_1(self):
        s = 'abb'
        #    -00
        z = z_array(s)
        self.assertEqual([3, 0, 0], z)

    def test_z_2(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_z_3(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_n_1(self):
        s = 'abb'
        #    01-
        n = n_array(s)
        self.assertEqual([0, 1, 3], n)

    def test_n_2(self):
        s = 'abracadabra'
        #    1004010100-
        n = n_array(s)
        self.assertEqual([1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11], n)

    def test_n_3(self):
        s = 'abababab'
        #    0204060-
        n = n_array(s)
        self.assertEqual([0, 2, 0, 4, 0, 6, 0, 8], n)

    def test_big_l_prime_1(self):
        s = 'abb'
        #    001
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 2], big_l_prime)

    def test_big_l_prime_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)

    def test_small_l_prime_1(self):
        s = 'abracadabra'
        # N  1004010100-
        # l'           1
        # l'        4
        # l' 44444444111
        small_l_prime = small_l_prime_array(n_array(s))
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

    def test_good_suffix_match_mismatch_1(self):
        p = 'GGTAGGT'
        big_l_prime, big_l, small_l_prime = good_suffix_table(p)
        self.assertEqual([0, 0, 0, 0, 3, 0, 0], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 3, 3, 3], big_l)
        self.assertEqual([7, 3, 3, 3, 3, 0, 0], small_l_prime)
        self.assertEqual(0, good_suffix_mismatch(6, big_l_prime, small_l_prime))
        self.assertEqual(0, good_suffix_mismatch(6, big_l, small_l_prime))
        #  t:      xT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(5, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(5, big_l, small_l_prime))
        #  t:     xGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(4, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(4, big_l, small_l_prime))
        #  t:    xGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(3, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(3, big_l, small_l_prime))
        #  t:   xAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(2, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(2, big_l, small_l_prime))
        #  t:  xTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(1, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(1, big_l, small_l_prime))
        #  t: xGTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(0, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(0, big_l, small_l_prime))

    def test_good_suffix_table_1(self):
        s = 'abb'
        #    001
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 2], big_l_prime)
        self.assertEqual([0, 0, 2], big_l)
        self.assertEqual([3, 0, 0], small_l_prime)

    def test_good_suffix_table_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        # l' -4444444111
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8], big_l)
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)


# My Enhanced Version of Python code
def ReadFasta(filename):
    with open(filename, 'rt') as input:
        t = ''
        for line in input:
            if '>' in line:
                continue
            else:
                t += line.strip('\n')
    return t
t = ReadFasta('chr1.GRCh38.excerpt.fasta')
# print(len(t))

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


def naive_with_counts(p, t):
    '''this function is used to naive mathcing pattern in text'''
    occurrences = []
    num_alignments = 0
    num_character_comparisons  = 0

    for i in range(len(t) - len(p) + 1):
        num_alignments += 1
        match = True # Set original flag
        for j in range(len(p)):
            if not t[i + j] == p[j]:
                match = False
                break
        
        if match:
            occurrences.append(i)
    num_character_comparisons = num_alignments * len(p)
    return occurrences, num_alignments, num_character_comparisons


def boyer_moore_with_counts(p, p_bm, t):
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(p) - 1, -1, -1):
            num_character_comparisons += 1
            if not p[j] == t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_alignments, num_character_comparisons




'''Q1: How many alignments does the naive exact matching algorithm try when matching 
the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to 
the excerpt of human chromosome 1?  (Don't consider reverse complements.)'''
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# print(len(p))
# print(len(t))
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(f'The naive exact matching try {num_alignments} times of alignment')
# 56922


'''Q2: How many character comparisons does the naive exact matching algorithm try 
when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to 
the excerpt of human chromosome 1?  (Don't consider reverse complements.)'''
print(f'The naive exact matching does {num_character_comparisons} times of character comparision')


'''Q3: How many alignments does Boyer-Moore try 
when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to 
the excerpt of human chromosome 1?  (Don't consider reverse complements.)'''
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# print(t)
# print(len(t))
p_bm = BoyerMoore(p)
# print(p_bm)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences)
print(f'The Boyer_moore matching try {num_alignments} times of alignment')
print(f'The Boyer_moore matching try {num_character_comparisons} times of character comparison')


'''Q4: Index-assisted approximate matching. In practicals, we built a Python class called Index
Implement the pigeonhole principle using Index to find exact matches for the partitions. 
Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches (substitutions). 
We will use an 8-mer index.
'''
# import bisect


# class Index(object):
#     """ Holds a substring index for a text T """

#     def __init__(self, t, k):
#         """ Create index from all substrings of t of length k """
#         self.k = k  # k-mer length (k)
#         self.index = []
#         for i in range(len(t) - k + 1):  # for each k-mer
#             self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
#         self.index.sort()  # alphabetize by k-mer, this is going to take a long time for a long text t

#     def query(self, p):
#         """ Return index hits for first k-mer of p """
#         kmer = p[:self.k]  # query with first k-mer
#         i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
#         hits = []
#         while i < len(self.index):  # collect matching index entries
#             if self.index[i][0] != kmer:
#                 break
#             hits.append(self.index[i][1])
#             i += 1
#         return hits


# def approximate_match(p, t, n):
    
#     segment_length = round(len(p) // (n+1))
#     all_matches = set()
#     for i in range(n+1):
#         start = i*segment_length
#         # print(f'The start is {start}')
#         end = min((i+1)*segment_length, len(p))
#         p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
#         matches, _, _ = boyer_moore_with_counts(p[start:end], p_bm, t)
#         # print(matches)

#         for m in matches:
#             if m < start or m-start+len(p) > len(t): # not to let p run off at the beginning or the end of t
#                 continue
            
#             # Count the mismatches
#             mismatches = 0
#             for j in range(0, start):
#                 # print(f'The j is {j}')
#                 if not p[j] == t[m-start+j]:
#                     mismatches += 1
#                     if mismatches > n:
#                         break
#             for j in range(end, len(p)):
#                 if not p[j] == t[m-start+j]:
#                     mismatches += 1
#                     if mismatches > n:
#                         break
#             if mismatches <= n:
#                 all_matches.add(m - start)
#     return list(all_matches)

# # print(approximate_match(p, t, 2))
# p = 'GGCGCGGTGGCTCACGCCTGTAAT'

# # Let's assume that mismatch is not allowed in the seed sequence.
# def queryIndex(p, t, index):
#     k = index.k
#     occurrences = []
#     for i in index.query(p):
#         mismatch = 0  # being modified
#         for j in range(k, len(p)):
#             
#             if p[j] != t[i+j]:
#                 mismatch += 1
#                 if mismatch >= 2:
#                     break
#         occurrences.append(i)
#     return occurrences

# index = Index(t, 8)
# occurrences = queryIndex(p, t, index)
# print('Allowed 2 mismatches (not including seed length), there are {}'.format(', '.join(map(str, occurrences)).rstrip(',')))
# # print(p)
# # print(t[57056:57056+24])

# # print('Allowed 2 mismatches (not including seed length), there are {} alignments'.format(len(occurrences)))


# Let's assume that mismatch is allowed in the seed sequence.
# 1. modify the Index class first
import bisect


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            # print(i)
            mismatch = 0
            for j in range(self.k):
                if p[j] != self.index[i][0][j]:
                    mismatch += 1
                if mismatch > 2:
                    break
            if mismatch <= 2:
                tmp = [self.index[i][1], mismatch]
                hits.append(tmp)
            i += 1
        return hits
        #     if self.index[i][0] != kmer:
        #         break
        #     hits.append(self.index[i][1])
        #     i += 1
        # return hits

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
index = Index(t, 8)
hits = index.query(p)
# print(hits[:5])
# print(p)
# print(t[662865:662873])
# print(t[542009:542009+8])
# print(t[57056:57056+24])

def queryIndex(p, t, hits):
    k = index.k
    occurrences = []
    for item in hits:
        mismatch = item[1]  # retrieve the original mismatches number
        # # 
        # if item[0] == 57056:
        #     print(item[1])
        # # 
        for j in range(k, len(p)):
            # # 
            # print(k)
            # #
            if p[j] != t[item[0]+j]:
                mismatch += 1
                if mismatch > 2:
                    break
        # print(mismatch)
        if mismatch <= 2:
            occurrences.append(item[0])
    return occurrences

occurrences = queryIndex(p, t, hits)
# print(occurrences[:5])
print('Allowed 2 mismatches (including seed length), there are {} alignments'.format(len(occurrences)))
# print(p)
# print(t[56922:56922+24])
# print(t[57056:57056+24]) # Bug results
# print(t[84641:84641+24])


'''Q5: Using the instructions given in Question 4, 
how many total index hits are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?
Note: My resutls are already in Q4
'''

'''Q6:
1) Background
Let's examine whether there is a benefit to using an index built using subsequences of T rather than substrings, as we discussed in the "Variations on k-mer indexes" video.  
We'll consider subsequences involving every N characters. For example, if we split ATATAT into two substring partitions, we would get partitions ATA (the first half) and TAT (second half).  
But if we split ATATAT into two  subsequences  by taking every other character, we would get AAA (first, third and fifth characters) and TTT (second, fourth and sixth).

Another way to visualize this is using numbers to show how each character of P is allocated to a partition.  Splitting a length-6 pattern into two substrings could be represented as 111222, 
and splitting into two subsequences of every other character could be represented as 121212
The following class SubseqIndex is a more general implementation of Index that additionally handles subsequences. It only considers subsequences that take every Nth character:

2) Question
Write a function that, given a length-24 pattern P and given a SubseqIndex object built with k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches.
When using this function, how many total index hits are there when searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?  
(Again, don't consider reverse complements.)
'''
import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)   # the total length of kmer is n, and from the start pos (assume 1), it is going to be (k-1) * section_length to extend
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

index = SubseqIndex(t, 8, 3)

# Let's assume there is no mismatches in alignment of seed sequence
def queryIndex(p, t, index):
    occurrences = []
    for i in index.query(p):
        mismatch = 0

        for j in range(len(p)):
            if p[j] != t[i+j]:
                mismatch += 1
            if mismatch > 2:
                break
        if mismatch <= 2:
            occurrences.append(i)
    # print(index.query(p))
    # print(occurrences)
    return occurrences

p = 'GGCGCGGTGGCTCACGCCTGTAAT'
occurrences = queryIndex(p, t, index)
# print(occurrences[:5])
print(p)
# # print(t[16614:16614+len(p)]) # Bug result
# # print(t[56922:56922+len(p)])
print(t[84641:84641+len(p)])


# Let's assume there is no mismatches in alignment of seed sequence
# Skip