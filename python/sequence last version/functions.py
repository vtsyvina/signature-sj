import timeit
from os import listdir
from os.path import isfile, join

from coursera import Profile


def sequence_coordinates(sequence):
    result = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for letter in sequence:
        result[letter] += 1
    return result

#checks if profile contains 1 of i's position
def is_equal_letter(profile, i):
    for letter in profile:
       if profile[letter][i] == 1:
           return True
    return False

#if profile contains more than one 1.000 value followed by each other,
# this functions will replace every equal substring with one letter
def reduce_letter_count(seq):
    sequences = [x for (_, x) in seq.items()]
    profile = Profile(sequences)
    for l in profile:
        profile[l].insert(0,0)
    start = -1
    for i in reversed(range(0,len(profile['A']))):
        if start == -1 and is_equal_letter(profile, i):
            start = i
        if not is_equal_letter(profile, i):
            if start != -1 and start-i > 1:
                for name in seq:
                    # we add 0 at 0 position, so no +1 in this case
                    seq[name] = seq[name][:i+1]+seq[name][start:]
            start = -1
    return seq


def measure_methods_time(k_range, files_count, mult, mut):
    folder = "cleaned_independent_264"
    onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]
    attempts = 1
    algs = ['points_method']#, 'brute_force_close_pairs']

    files = ["VNHCV41_1b_"+str(mult)+"_"+str(mut)+".fas",
             "VNHCV55_1b_"+str(mult)+"_"+str(mut)+".fas",
             "VNHCV227_1b_"+str(mult)+"_"+str(mut)+".fas"]
    test_files= ['AMC_P08_1a.fas']
    for file in test_files:#onlyfiles[:files_count]:
        for k in k_range:
            for algorithm in algs:

                statement = algorithm+'(sequences, k)'
                setup = '''from close_pairs_algorithms import {}
from utils import read_sequences_with_additional_info
from utils import read_sequences
folder = "cleaned_independent_264/"
sequences = read_sequences(folder+\"{}\")
l = len(sequences)
print("all pairs "+str(l*(l-1)//2))
k = {}'''.format(algorithm, file, k)
                # if algorithm != 'brute_force_close_pairs' or (k == 3 and 'VNHCV227_1b' in file):
                print("file: {:<15} alg: {:<40} k: {} time: {}".format(file, algorithm, k, str(timeit.timeit(statement, number=attempts, setup=setup) / attempts)))
            print()

# process input files if they are not in proper format
def clean_data():
    folder = "cleaned_independent_264"
    onlyfiles = [folder+"/"+f for f in listdir(folder) if isfile(join(folder, f))]
    for filename in onlyfiles:
        res = []
        with open(filename, "r") as input_file:
            name = ''
            seq = ''
            for line in input_file:
                tmp = line.rstrip()
                if tmp.startswith(">") :
                    name = tmp
                elif len(tmp) == 0:
                    res.append((name, seq))
                    seq = ''
                else:
                    seq+= tmp
        output = open(filename, "w")
        for (name, seq) in res:
            output.write(name+"\n")
            output.write(seq+"\n")

# builds points dict by given sequences
def build_points(sequences):
    points = {}
    for name, sequence in sequences.items():
        coords = sequence_coordinates(sequence)
        coords = (coords['A'], coords['C'], coords['G'], coords['T'])
        if coords not in points:
            points[coords] = []
        points[coords].append(
            (name, sequence, coords))
    return points

def build_sequences_kmer_hashes(sequences, length):
    kmers = {}
    l_to_d = {'A' : 0, 'C': 1, 'G': 2, 'T': 3}
    exp = 4 ** (length-1)
    for name, sequence in sequences.items():
        kmers[name] = {}
        hash_value = 0
        for letter in sequence[:length]:
            hash_value *= 4
            hash_value += l_to_d[letter]
        kmers[name][hash_value] = 1
        kmers[name]["coherently0"] = hash_value
        for i in range(1, len(sequence)-length+1):
            hash_value -= exp * l_to_d[sequence[i-1]]
            hash_value *= 4
            hash_value += l_to_d[sequence[i+length-1]]
            kmers[name][hash_value] = 1
            if i % length == 0:
                kmers[name]["coherently"+str(i//length)] = hash_value
    return kmers

# Builds dictionary with folloving form:
#   123213213 -> set(seq_1_name, seq_2_name, ..., seq_17_name) - set of all sequences in sample that contain string with given hash
#   342342344 -> set(seq_7_name, seq_3_name, ..., seq_10_name)
#   seq_3_name -> [3232312,453412384, 321315432, ...] - list of string hashes in sequence on 0, 1*lenght, 2*length.... places
#   whole_sample -> [set(3232312,453412384, 321315432), set(9937312,453492384, 956315432)] - sets of all posible hashes for given position in sample
def build_kmers_dict(sequences, length):
    kmers = {}
    l_to_d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    exp = 4 ** (length - 1)
    kmers["whole_sample"] = [set() for i in range(len(list(sequences.values())[0])//length)]
    for name, sequence in sequences.items():
        hash_value = 0
        for letter in sequence[:length]:
            hash_value *= 4
            hash_value += l_to_d[letter]
        if hash_value not in kmers:
            kmers[hash_value] = set()
        kmers[hash_value].add(name)
        kmers[name] = []
        kmers[name].append(hash_value)
        kmers["whole_sample"][0].add(hash_value)
        for i in range(1, len(sequence) - length + 1):
            hash_value -= exp * l_to_d[sequence[i - 1]]
            hash_value *= 4
            hash_value += l_to_d[sequence[i + length - 1]]
            if hash_value not in kmers:
                kmers[hash_value] = set()
            kmers[hash_value].add(name)
            if i % length == 0:
                kmers[name].append(hash_value)
                kmers["whole_sample"][i // length].add(hash_value)
    return kmers

def lower_bound_estimate(coord1, coord2):
    positive = 0
    negative = 0

    for letter in 'ACTG':
        dif = coord1[letter] - coord2[letter]
        if dif > 0:
            positive += dif
        else:
            negative += dif
    return positive if positive > -negative else -negative


def lower_bound_estimate_list(coord1, coord2):
    positive = 0
    negative = 0
    for coord in range(len(coord1)):
        dif = coord1[coord] - coord2[coord]
        if dif > 0:
            positive += dif
        else:
            negative += dif
    return positive if positive > -negative else -negative

# takes dictionary of sequences name -> sequence, step(piece positions), piece lenght
# returns dict of all possible pieces with analogical dictionary as input
def build_segments_dict(sequences, step, piece_lenght):
    result = {}
    for (name, seq) in sequences.items():
        segment = seq[step*piece_lenght:(step+1)*piece_lenght];
        if segment not in result:
            result[segment] = {}
        result[segment][name] =  seq
    return result

def kmer_validation(hashes1, hashes2, coincidences, hashes_count):
    count = 0
    for i in range(hashes_count):
        if hashes1["coherently"+str(i)] in hashes2:
            count += 1
        if hashes_count - i < coincidences - count:
            return False
    return True
