# text = 'AGTAGCCACTAGCCACTAGCCACTAGCCACTCCCATGAGAGCCACTCAGCCACTAGCCACTTAGCCACTTAGCCACTCCCGGCCAAGCCACTGAGCCACTGTCAATAGCCACTTCCTAGCCACTAGCCACTCTCAAGCCACTAAGCCACTAGTAGCCACTTACAAGCCACTAGCCACTAGCCACTAGCGAGCCACTAGCCACTAGCCACTTAAAGCCACTCCCGAGCCACTTCCATTAGCCACTAGCCACTAAGCCACTTTCAGCCACTAGCCACTCACACAGCCACTACAAAAATGTATGAGCCACTTAGAAGCCACTAGCCACTATCTAGCCACTTGAGCCACTTCCCAGCCACTTGAGCCACTAAAGCCACTTAGCCACTTACGTGAGCCACTCCAGCCACTGAGCCACTGAAGCCACTAAAGCCACTCATATAGCCACTTAGCCACTTACAGCCACTAGCCACTCAGAGCCACTAGCCACTGAGCCACTGGTAGCCACTAGCCACTTGCATTGCAGCCACTCAGCCACTGAGCCACTGGGTGATCTCCAGCCACTCGTCAGTAGCCACTAGAGCCACTAAGTATGCATCAGCCACTAGCCACTGTAGAGCCACTATAGTTCTCCGCTGGACCAGCCACTTAGCCACTGTTCGCTCGTTAGCCACTAGCCACTTAGCCACTAGCCACTGAGCCACTGCCCGGATAACCCAGCCACTGTAGCCACTCTTACGAGCCACTCCATGGCCATGGCAGCCACTGAGCCACTTTGATTAGCCACTAAAGCCACTCCAGCCACTATGAGCCACTGGAAGGAGCCACTCAAGCCACTGTAGCCACTAGCCACTGTCCGCATTAGAAAAAAATCTGTCAGCCACTTTCTGGTAGCCACTAGCCACTTTAAACAGCCACTTAGCCACTGAAAGAGCCACTGAGCCACTAAAGCCACTTCGAGCCACTCGATAGCCACTGAGCCACTAGCCACTAGCCACTA'
# pattern = 'AGCCACTAG'

# def PatternCount(text,pattern):
#     count = 0
#     for i in range(len(text)-len(pattern)+1):
#         if text[i:i+len(pattern)]==pattern:
#             count = count+1
#     return count
# print(PatternCount(text,pattern))

# def frequentpattern (seq,k):
#     frequency_map = {}           #initialising empty dictionary
#     n = len(seq)
#     for i in range (n-k+1):     # This creates a loop that will repeat the following code for the number of times indicated between parenthesis
#         pattern = seq[i:i+k]           #code4-10 = fiiling dictionary 
#         if pattern not in frequency_map :
#              frequency_map[pattern]=1
#         else:
#              frequency_map[pattern] +=1
#     return frequency_map
# e = frequentpattern("CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT",3)  
# print(e)
# maximum_value = max(e.values())                       #finding maximum value using max function 
# print(maximum_value)
# the_most_frequent_pattern = [x for x in e if e[x]==maximum_value]       #list comprehension 
# print(the_most_frequent_pattern)

# from Bio.Seq import Seq
# sequence = Seq('')
# reverse_complement = sequence.reverse_complement()
# print (reverse_complement)
# genome = 'ATGACTTCGCTGTTACGCGC'
# pattern = 'CGC'
# position = []
# for i in range (len(genome)-len(pattern)+1):
#     if genome[i:i+len(pattern)]==pattern:
#         position.append(i)
# print (position)

# import numpy

# def FindClumps(text, k, L, t):
#     n = len(text)
#     out = []
#     # FOR each window of size L
#     for i in range(n-L):
#         pattern = []
#         window = text[i:i+L]
#         # FOR each pattern present in window
#         for j in range(len(window)-k):
#             pattern.append(window[j:j+k])
#         # find unique k-mers and their frequencies
#         kmers, counts = numpy.unique(pattern, return_counts=True)        
#         out += kmers[counts>=t].tolist()
#     return numpy.unique(out, return_counts=True)

# with open('ecoli.txt', 'r') as f:
#     text = f.read()
# f.close()

# print (FindClumps(text,9,29,4))


# def hummingDistance(text1,text2):
#     count = 0
#     for i in range(len(text1)):
#         if text1[i] == text2[i]:
#             count == count 
#         else:
#             count = count + 1
#     return count
# print (hummingDistance('TAGCCCTAGTGGTTTGCGAACACTAGAACCAGTAAAGCAAAAAAAACCCCCCTTGCGCAGTCATAAGTCGCTTTTAGCCCCTCGCACAGTTCTGAGTATTCACGGATATTGAGTCTTCAGAAAGGAGTCATGTCCTGGTTCAGTGGTCGGCTACGGGGTTGCAGGTCAGCATTATAGTCATCAACTAATGGTCGCCATTCCGGAACGGAAGATACACTTGGAAGCTACGGTCGGGTGAGACCAAGCTCACCAATCTAGCGCCCTAATAGTTTGTCTTTAAAATGCTTTTTGTGTTACAGAAGCTAGTTCGTTCCCGGTCCTCGACAGCCGGTGTTTTCTTTAGCACCTTTCTTCAGATGAACGGTCCGCCTCTGGTACGCCCTCCACTAGTATGCTTGCATTCTCGCCTTTAGGGTGCCGCTGTACCTAGACTGCAGAGCGCGAAATTGTATCATTCCATAACGAGGCTAATTGGTCATGGATAACTCTTATGATGGGCCGGGGAGGGAGTCCCTTTGAGACCAGTTGCAGGAAAGACTTTCTACAGTAAGCCTGGTTAGACGCCCTACCAGATACAGGTGTGTAATCCGTATTCATTACGTGAGGATTGTCTTGGGGGTTAATCGATACAATATTTCACATTAGACGCGGTAACGTGTGTAACTTCGTAGCCGTAGCGCTTCTCCACGTATTGCTAAGTTCGGTCTATATCCACGTACTGATAAATCAGTTTCTGTAACGTGCTCCCCTGCTTTTAAACACGCTTTGTATCTTGTAATGAACGAGATGCGTCCCAGATCAGTTGGGCTTGAGAGAAGGCATCTTGCCGTAACGACAAGTCCGACATTGCAGGCAACGTACTGCGTCACAAGGCCTCTGTTCTGGAGAATGCTTTCAGTCTTGAGCTCAGGAAGAGTCTACGTTCTCTTGTCCCGCTCAACGACCGTCGGGAAGGCCTGACGTCGACGATTTAAACTGACGGGCAAGCATTTTTGCAAGATGCTACTGAGCATACTACGTCTGGTTGGTTCCGAATCTTAGCACCGATTAGCAGCGCATA','TCGTTACCCAACAAACCGGGATTACGAGGGTGATGCCAAGACTTTTGGTGCGGTCTTCCCCGCCAGTGGACAGGCGTCGGGCTATCTGATGTTAAATACTTTATATGCAGACAGGGGCACGTCATCGGGGCCGCGTCAAGTATTAATCTTCACCGCGAGATCCCAACCCTTAAAACTTCCTGACCTACCGAACTCCTTCCGTGGCGCAGACCCGAGGAGTCCCCGGCTAATCACACCCTAGTCATAGACCGACTTAGCCTCAACAACTAAAAAGAAAGTCTACACTAACTTAATGATGACCAAAGGACTAGTATCGTTTGACTAATTGGTAATGCAACGAGGTTCACTCTCATATACCTTTTGATCACATCCGGCACGTATTTATTGATTGAAATAGGTAGGATGACGCCCGAATCCGACGCTAGGGGCGGCGGGCTCTTCGTGACATACTTATCGGCGATTCTGCATTTATGTATGCGTTTCACAAGGGGCAGCTAGTATCCGGACGTTATGATTTAAGCCTATGTGGCTCTAAGTGCTGCAATTAGCGTTTTCCGGCGGCGAAGCGAATTGCCACTCCACCAGGATCTGAGCTGGAGTACACATCTCATTTGGAGGAGGGTTACCCCTTGAAGGGCCTCGGCTCATCACTAGCGTACGACATATAGCAGTAACGAGCCGTTCATATTGGTAAAGCTCTGCCAGGCGAGACACGGGTGGTCCAGTCCACGAAATCTCCTGGAACCAGTTTCTCGGTTAGCAGACCGGAGGACTCCCCCCGGCCCCTGAATAGTAGATGGAAGCCGTTACACGGCCACCGTTAAATCGCTACTAATGTTTATGCCAAGCATAATGCACTACGCGGCCTATCAGCGAGCCGTGGAATCGGTACCTCTGCACCTAGTAGCCCGATCCTTGAGTGGATTGACGAAATAACGCCAGCCACAGAGACTGGACGGAGAAACGGGTGGGGTTGTTCCTACTTGTCTCCCTGAGGTCGACTTTTTTTTGCACTGCCTGTTGAACAGTTTACGTTTGACTAACTGAGGCGTTACCTCTT'))

# with open('dataset_7_10.txt', 'r') as f:
#     text = f.read()
# f.close()

# def skew(text):

#     skew = 0

#     skew_list= []

#     skew_list.append(0)

#     for i in range(len(text)):

#         if text[i] == 'C' :

#             skew = skew - 1

#         elif text[i] == 'G':

#             skew += 1

#         else :

#             skew = skew

#         skew_list.append(skew)

#     return skew_list
# print (skew(text))

# def min_skew(text):

#     skew_list = skew(text)

#     minimum = min(skew_list)

#     min_skew_list = list()

#     for i in range(len(text)):

#         if skew_list[i] == minimum:

#             min_skew_list.append(i)

#     min_skew_list = str(min_skew_list)

#     min_skew_list = min_skew_list.replace(',','')

#     print(min_skew_list)


# print(min_skew(text))


# position = []
# def approx_pat(pattern, genome, d):
#     for i in range (len(genome) - len(pattern)+1):
#         if hamming_distance(pattern, genome[i:i + len(pattern)]) <= int(d):
#             position.append(i)
#     return position
# def hamming_distance(q, p):
#         dist = 0
#         for i in range(len(p)):
#             if p[i] != q[i]:
#                 dist += 1
#         return dist
# print(len(approx_pat('ACCGGTT', 'TGCTCTTGAGGTTGAGTGCGACCGGTTATCGACCTCGCAACACCTCGTGCGAATAAAGCCACACTTCTTTTGACGATTTGATCGTCATGAGACTGGACCAAAAAACAGTGTGGCTTGAGGTTGATTAAAACTTAAGATGTGAGGGGGTCCCAGAAGGACCCCAACTCCCGCTTCCTTAGGCTCCACTGACTCACCTGCCAATGAAGTGGCTTCAGGGGACTATTCCGGGTTGTGTGATCATCCGCCAAGAACATTGATGCGACACTGGTAATAAATAGGAGACGATCCGCCGGTTGCAGCCAACGGATATCTTCCGACTACTAAGATGGACAAGATGTGCAGAAGAG', '3')))


# import itertools
# import time
# from collections import defaultdict

# def FrequentWordsWithMismatches(Genome, k, d):
#     start = time.process_time()
#     aprox_frq_words = []
#     frequencies = defaultdict(lambda: 0)
#     # all existent kmers with d mismatches of current kmer in genome
#     for index in range(len(Genome) - k + 1):
#         curr_kmer_and_neighbors = PermuteMotifDistanceTimes(Genome[index : index + k], d)
#         for kmer in curr_kmer_and_neighbors:
#             frequencies[kmer] += 1 

#     for kmer in frequencies:
#         if frequencies[kmer] == max(frequencies.values()):
#             aprox_frq_words.append(kmer)
#     end = time.process_time()
#     print("Time:", end - start)
#     return aprox_frq_words


# def PermuteMotifOnce(motif, alphabet={"A", "C", "G", "T"}):
#     """
#     Gets all strings within hamming distance 1 of motif and returns it as a
#     list.
#     """

#     return list(set(itertools.chain.from_iterable([[
#         motif[:pos] + nucleotide + motif[pos + 1:] for
#         nucleotide in alphabet] for
#         pos in range(len(motif))])))


# def PermuteMotifDistanceTimes(motif, d):
#     workingSet = {motif}
#     for _ in range(d):
#         workingSet = set(itertools.chain.from_iterable(map(PermuteMotifOnce, workingSet)))
#     return list(workingSet)

# print (FrequentWordsWithMismatches('GGCTGAACGGCTCGCGCGGTTTGCGGAACGCGGTTTGGCTCGCGTTTGAACGAACGGCTGGCTGTTTGTTTGCGGAACGGCTGAACGTTTGCGGTTTCGCGTTTCGCGCGGAACGTTTGCGGGCTGCGGAACGAACGTTTGGCTGGCTGTTTGCGCGCGAACGAACCGCGAACGAACGTTTGTTTGCGGTTTCGCGGCTGGCTGCGGAACGCGGAACGTTTCGCGCGGCGGCGCGCGAACGGCTGGCTGGCTGTTTCGCGAACGTTTGTTTGAACCGCGCGGGCTGAACGAACGTTTGGCTGCGCGCGCGGGCTGAACGTTTCGC', 5, 3))




# def mismatches_and_rc(text,k,d):
#     frequencyarray={}
#     rev_text=rev_comp(text)
#     for i in range(0,4**k):
#         frequencyarray[i]=0
#     for i in frequencyarray:
#         count,r_count=0,0
#         count=PatternMatching(text,NumberToPattern(i,k),d)
#         r_count=PatternMatching(rev_text,NumberToPattern(i,k),d)
#         frequencyarray[i]=count+r_count
#     return frequencyarray
            
        
# def PatternMatching(text,pat,d):
#     count=0
#     for i in range(0,len(text)):
#         p=text[i:i+len(pat)]
#         if len(p)!= len(pat): break
#         if hammingdistance(p,pat)<=d:
#             count=count+1
#     return count

# def PatternToNumber(pattern):
#     if len(pattern)==0: return
#     SymbolToNumber = {'A':0,'C':1,'G':2,'T':3}
#     if len(pattern)==1: return SymbolToNumber[pattern]
#     n=len(pattern)
#     symbol=pattern[n-1]
#     prefix=pattern[:n-1]
#     return (4*PatternToNumber(prefix)+SymbolToNumber[symbol])

# def NumberToPattern(index,k):
#     NumberToSymbol = {0:'A',1:'C',2:'G',3:'T'}
#     if k==1: return NumberToSymbol[index]
#     prefix_index=index//4
#     r=index%4
#     symbol=NumberToSymbol[r]
#     return NumberToPattern(prefix_index,k-1)+symbol

def hammingdistance(p,q):
    count=0
    if len(p)!= len(q):
        print("The two sequences vary in lenghts")
        return
    for i in range(0,len(p)):
        if p[i]!=q[i]:
            count=count+1
    return count
print (hammingdistance('CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG','ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'))

# def rev_comp(pattern):
#     rpattern=''
#     a=[]
#     for i in pattern:
#         if i=='A':
#             a.append('T')
#         if i=='T':
#             a.append('A')
#         if i=='C':
#             a.append('G')
#         if i=='G':
#             a.append('C')
#     rpattern =rpattern.join(a)
#     return(rpattern[::-1])


# text=input('Enter the text')
# k=int(input('Enter length of k-mer'))
# d=int(input('Enter d'))
# approx_kmers= mismatches_and_rc(text,k,d)
# for a in approx_kmers:
#     if approx_kmers[a]== max(approx_kmers.values()):
#         print(NumberToPattern(a,k),end=' ')

# def immediateNeighbors(pattern):
#     neighborhood = set()
#     for i in range (len(pattern)):
#         symbol = pattern[i]
#         for x in 'ACTG':
#             if x != symbol:
#                 neighbor = pattern[:i] + x + pattern[i+1:]
#                 neighborhood.add(neighbor)
#     return neighborhood
# print (immediateNeighbors())



chars = "ACGT"
def neighbors(pattern, d):
    assert(d <= len(pattern))
    if d == 0:
        return [pattern]
    r2 = neighbors(pattern[1:], d-1)
    r = [c + r3 for r3 in r2 for c in chars if c != pattern[0]]
    if (d < len(pattern)):
        r2 = neighbors(pattern[1:], d)
        r += [pattern[0] + r3 for r3 in r2]
    return r

print (len(neighbors("ACGT", 3)))