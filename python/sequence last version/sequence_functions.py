import itertools
import timeit
from os import listdir
from os.path import isfile, join
from time import time

import Levenshtein

from close_pairs_algorithms import brute_force_close_pairs, equal_kmer_method, equal_kmer_method_two_samples, \
    brute_force_close_pairs_two_samples
from close_pairs_algorithms import points_method
from close_pairs_algorithms import points_method_two_samples, points_method_with_kmer_hashes
from coursera import Profile, Consensus
from functions import reduce_letter_count, build_points, \
    build_sequences_kmer_hashes, build_kmers_dict
from utils import read_sequences

# count = 10
# ed = timeit.timeit("sequence_coordinates(s)", number=count, setup='''from functions import sequence_coordinates
# input_file = open("VNHCV227_1b.fas","r")
# input_file.readline()
# s = input_file.readline().rstrip()
#
# input_file.readline()
# t = input_file.readline().rstrip()
# ''')
# print (str(ed/count))
#
# ed = timeit.timeit('''result = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
# for letter in s:
#     result[letter] += 1''', number=count, setup='''from functions import sequence_coordinates1
# input_file = open("VNHCV227_1b.fas","r")
# input_file.readline()
# s = input_file.readline().rstrip()
#
# input_file.readline()
# t = input_file.readline().rstrip()
# ''')
# print (str(ed/count))
#
count = 100
ed = timeit.timeit("distance(s, t)", number=count, setup='''from Levenshtein import distance
input_file = open("VNHCV227_1b.fas","r")
input_file.readline()
s = 'CACCGACTGCGGCGCTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGAACATCGCGCTTATCA'

input_file.readline()
t = 'CACCGACTGCGGCACTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGGACATCGCGCTTATCA'
''')
print (str(ed/count))
#
# ed = timeit.timeit("hamming(s, t)", number=count, setup='''from Levenshtein import hamming
# input_file = open("VNHCV227_1b.fas","r")
# input_file.readline()
# s = input_file.readline().rstrip()
#
# input_file.readline()
# t = input_file.readline().rstrip()
# ''')
# print (str(ed/count))
def test():
    count = 0
    folder = "cleaned_independent_264"
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]
    d = {}
    r = {}
    r['c'] = "ACTTAGAAC"
    r['d'] = "AGTTACCAC"
    e = reduce_letter_count(r)
    allseq = []
    for file in onlyfiles:
        seq = read_sequences("cleaned_independent_264/" + file)
        sequences = [x for (_, x) in seq.items()]
        if len(sequences[0]) == 264:
            allseq.extend(sequences)
        # print_profile(Profile(sequences))
        start_len = len(sequences[0])
        seq2 = reduce_letter_count(seq)
        sequences = [x for (_, x) in seq.items()]
        print("TOTAL PROFILE")

        if (len(sequences) > 200):
            profile = Profile(sequences)
            count = 0
            all = 0
            prev = False

            for i in range(len(profile['A'])):
                cur = False
                for letter in 'ACTG':
                    if profile[letter][i] == 1:
                        all += 1
                        cur = True
                if prev and cur:
                    count +=1
                prev = cur
            print("{:<20} {:<3} {:<3} {:<3} {:.3f}".format(file, all,len(profile['A']),start_len, len(profile['A'])/start_len))
    print_profile(Profile(allseq))


def test2():
    folder = "cleaned_independent_264"
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]
    con = []
    start = time()
    for file in onlyfiles:
        print(len(read_sequences("cleaned_independent_264/" + file)))
    print(time()-start)

def print_profile(profile):
    for letter in 'ACTG':
        print(letter+" ", end='')
        for s in profile[letter]:
            if s != 1:
                print(" {0:.2f} ".format(s), end='')
            else:
                print("*{0:.2f}*".format(s), end='')
        print()


def test3():
    folder = "cleaned_independent_264"
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]
    all_seq = {}
    points = {}
    hashes = {}
    dicts = {}
    whole_time_start = time()
    files_sub_set = onlyfiles
    # s = read_sequences(folder + "/AMC_P43_1a.fas")
    # t = read_sequences(folder + "/VAO_P48_1a.fas")
    # l = list(s.items())
    # left = l[:len(l)//2]
    # right = l[len(l)//2:]
    # s = {key: v for (key,v) in left}
    # t = {key: v for (key,v) in right}
    # p1 = build_points(s)
    # p2 = build_points(t)
    # hash_start = time()
    # points_method_two_samples(s,t,p1,p2,3)
    # print("POINTS METHOD!!!! " + str(time() - hash_start))
    # k1 = build_kmers_dict(s, 15)
    # k2 = build_kmers_dict(t, 15)
    # hash_start = time()
    # equal_kmer_method_two_samples(s, t, k1, k2, [], [], 3)
    # print("equal_kmer METHOD!!!! " + str(time() - hash_start))
    # hash_start = time()
    # brute_force_close_pairs_two_samples(s,t, 3)
    # print("brute_force_close_pairs_two_samples METHOD!!!! " + str(time() - hash_start))
    for f in files_sub_set:
        all_seq[f] = read_sequences(folder+"/"+f)
        print(f)
        # points_method(all_seq[f], 10)
        points[f] = build_points(all_seq[f])
        # hashes[f] = build_sequences_kmer_hashes(all_seq[f], 15)
        hash_start = time()
        dicts[f] = build_kmers_dict(all_seq[f], 14)
        print("Hash time " + str(time() - hash_start))
    found = []
    c = 0
    f1 = 'AMC_P01_1b.fas'
    f2 = 'AMC_P05_1b.fas'

    # close_pairs = equal_kmer_method_two_samples(all_seq[f1], all_seq[f2], dicts[f1], dicts[f2], points[f1], points[f2], 8)
    for (f1,f2) in itertools.combinations(files_sub_set, r=2):
        # if not (f1 == 'AMC_P04_1b.fas' and f2 == 'AMC_P05_1b.fas'):
        #     continue
        s = all_seq[f1]
        t = all_seq[f2]
        starttime = time()
        print("Start compairing "+f1+ " "+ f2)
        # close_pairs = points_method_two_samples(s, t, points[f1], points[f2], 8)

        close_pairs = equal_kmer_method_two_samples(s, t, dicts[f1], dicts[f2], points[f1], points[f2], 8)
        time1 = time() - starttime
        if time1 > 0.5:
            print("Time = " + str(time1)+" len "+str(len(close_pairs)))
        if len(close_pairs) > 0:
            print("!!!!!!!!! FOUND !!!!!")
            found.append((f1,f2))
        print()
    print("Whole time "+str(time()-whole_time_start))
    print("Related samples "+str(len(found)))
    for t in found:
        print(t)



# measure_methods_time(k_range=[3],files_count=70, mult=0, mut=0)
# generate_files()
test2()

