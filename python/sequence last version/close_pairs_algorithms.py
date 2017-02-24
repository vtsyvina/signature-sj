import itertools
from time import time

from Levenshtein import distance as lev_dis
from Levenshtein import hamming

from functions import sequence_coordinates, reduce_letter_count, lower_bound_estimate, lower_bound_estimate_list, \
    kmer_validation, build_kmers_dict
from utils import get_3_neighbourhood

debug = True


def brute_force_close_pairs(sequences, k):
    close_pairs = []
    count = 0
    for (p1, p2) in itertools.combinations(sequences.items(), 2):
        d = lev_dis(p1[1], p2[1])
        count += 1
        if d <= k:
            close_pairs.append((p1, p2, d))
    if debug:
        print(len(close_pairs))
        print("brute_force comps " + str(count))
    return close_pairs


def brute_force_close_pairs_two_samples(sequences1, sequences2, k):
    close_pairs = []
    count = 0
    for (p1, p2) in itertools.product(sequences1.items(), sequences2.items(), repeat=1):
        d = lev_dis(p1[1], p2[1])
        count += 1
        if d <= k:
            close_pairs.append((p1, p2, d))
    if debug:
        print("brute_force len " + str(len(close_pairs)))
        print("brute_force comps " + str(count))
    return close_pairs


def lower_bound_method(sequences, k):
    coordinates = []
    for name, sequence in sequences.items():
        coordinates.append((name, sequence, sequence_coordinates(sequence)))
    close_pairs = []
    comp = 0
    iters = 0
    for (p1, p2) in itertools.combinations(coordinates, 2):
        iters += 1
        if lower_bound_estimate(p1[2], p2[2]) <= k and hamming(p1[1], p2[1]) <= k * 2:
            d = lev_dis(p1[1], p2[1])
            comp += 1
            if d <= k:
                close_pairs.append((p1[0], p2[0], d))
    if debug:
        print("iters = " + str(iters))
        print("vomps = " + str(comp))
        print(len(close_pairs))
    return close_pairs


def points_method(sequences, k):
    points = {}
    reduce_letter_count(sequences)
    for name, sequence in sequences.items():
        coords = sequence_coordinates(sequence)
        coords = (coords['A'], coords['C'], coords['G'], coords['T'])
        if coords not in points:
            points[coords] = []
        points[coords].append(
            (name, sequence, coords))
    close_pairs = []
    iters = 0
    comp = 0
    to_process = set(points)
    reduce = 0
    diff_edit_hamming = 0
    for point in points:
        for coord in get_3_neighbourhood():
            sum_coord = coord[0] + coord[1] + coord[2]
            point2 = (point[0] + coord[0], point[1] + coord[1], point[2] + coord[2], point[3] - sum_coord)
            if point2 in to_process:
                iters += 1
                if point == point2:
                    for (p1, p2) in itertools.combinations(points[point], r=2):
                        iters += 1
                        h = hamming(p1[1], p2[1])
                        if h <= k:
                            close_pairs.append((p1[0], p2[0], h))
                        else:  # hamming(p1[1], p2[1]) <= 2 * k:
                            d = lev_dis(p1[1], p2[1])
                            comp += 1
                            if (d != h and h - d > diff_edit_hamming):
                                print("d " + str(d) + " h " + str(h))
                                diff_edit_hamming = h - d
                            if d <= k:
                                close_pairs.append((p1[0], p2[0], d))
                elif lower_bound_estimate_list(point, point2) <= k:
                    for (p1, p2) in itertools.product(points[point], points[point2], repeat=1):
                        iters += 1
                        h = hamming(p1[1], p2[1])
                        if h <= k:
                            close_pairs.append((p1[0], p2[0], h))
                        else:  # hamming(p1[1], p2[1]) <= 2 * k:
                            d = lev_dis(p1[1], p2[1])
                            comp += 1
                            if (d != h and h - d > diff_edit_hamming):
                                print("d " + str(d) + " h " + str(h))
                                diff_edit_hamming = h - d
                            if d <= k:
                                close_pairs.append((p1[0], p2[0], d))
        to_process.remove(point)
    if debug:
        print("points iters  = " + str(iters))
        print("points comps  = " + str(comp))
        print("reduce        = " + str(reduce))
        print("diff_edit_hamming= " + str(diff_edit_hamming))
        print(len(close_pairs))
    return close_pairs


def points_method_with_kmer_hashes(sequences, hashes, hashes_parts, k):
    points = {}
    reduce_letter_count(sequences)
    for name, sequence in sequences.items():
        coords = sequence_coordinates(sequence)
        coords = (coords['A'], coords['C'], coords['G'], coords['T'])
        if coords not in points:
            points[coords] = []
        points[coords].append(
            (name, sequence, coords))
    close_pairs = []
    iters = 0
    comp = 0
    to_process = set(points)
    reduce = 0
    diff_edit_hamming = 0
    for point in points:
        for point2 in to_process:
            iters += 1
            if point == point2:
                for (p1, p2) in itertools.combinations(points[point], r=2):
                    if kmer_validation(hashes[p1[0]], hashes[p2[0]], hashes_parts - k, hashes_parts):
                        iters += 1
                        h = hamming(p1[1], p2[1])
                        if h <= k:
                            close_pairs.append((p1[0], p2[0], h))
                        else:
                            d = lev_dis(p1[1], p2[1])
                            comp += 1
                            if d <= k:
                                close_pairs.append((p1[0], p2[0], d))
            elif lower_bound_estimate_list(point, point2) <= k:
                for p1 in points[point]:
                    for p2 in [p for p in points[point2] if
                               p1[0] != p[0] and kmer_validation(hashes[p1[0]], hashes[p[0]], hashes_parts - k,
                                                                 hashes_parts)]:
                        iters += 1
                        h = hamming(p1[1], p2[1])
                        if h <= k:
                            close_pairs.append((p1[0], p2[0], h))
                        else:
                            d = lev_dis(p1[1], p2[1])
                            comp += 1
                            if d <= k:
                                close_pairs.append((p1[0], p2[0], d))
        to_process.remove(point)
    if debug:
        print("points iters  = " + str(iters))
        print("points comps  = " + str(comp))
        print("reduce        = " + str(reduce))
        print("diff_edit_hamming= " + str(diff_edit_hamming))
        print(len(close_pairs))
    return close_pairs


def append_pair_if_needed(close_pairs, comp, iters, k, p1, p2):
    iters += 1
    h = hamming(p1[1], p2[1])
    if h <= k:
        close_pairs.append((p1[0], p2[0], h))
    else:
        d = lev_dis(p1[1], p2[1])
        comp += 1
        if d <= k:
            close_pairs.append((p1[0], p2[0], d))
    return comp, iters


def points_method_two_samples_hashes(sequences1, sequences2, points1, points2, hashes1, hashes2, hashes_parts, k):
    print("all pairs count = " + str(len(sequences1) * len(sequences2)))
    close_pairs = []
    iters = 0
    comp = 0
    reduced = 0
    throw = 0
    equal_length = len(list(sequences1.values())[0]) == len(list(sequences2.values())[0])
    numbers = {}
    for point1 in points1:
        for point2 in points2:
            iters += 1
            if lower_bound_estimate_list(point1, point2) <= k:
                for p1 in points1[point1]:
                    for p2 in [p for p in points2[point2] if
                               kmer_validation(hashes1[p1[0]], hashes2[p[0]], hashes_parts - k, hashes_parts)]:
                        iters += 1
                        if equal_length:
                            h = hamming(p1[1], p2[1])
                        else:
                            h = 1000000

                        if h <= k:
                            reduced += 1
                            close_pairs.append((p1[0], p2[0], h))
                        else:
                            d = lev_dis(p1[1], p2[1])
                            comp += 1
                            if d <= k:
                                close_pairs.append((p1[0], p2[0], d))
            else:
                throw += len(points1[point1]) * len(points2[point2])
    if debug:
        print("points iters  = " + str(iters))
        print("points comps  = " + str(comp))
        print("throw         = " + str(throw))
        print("numbers       = " + str(numbers))
        print(len(close_pairs))
    return close_pairs


def points_method_two_samples(sequences1, sequences2, points1, points2, k):
    print("all pairs count = " + str(len(sequences1) * len(sequences2)))
    close_pairs = []
    iters = 0
    comp = 0
    reduced = 0
    throw = 0
    equal_length = len(list(sequences1.values())[0]) == len(list(sequences2.values())[0])
    for point1 in points1:
        for point2 in points2:
            iters += 1
            if lower_bound_estimate_list(point1, point2) <= k:
                s = 2
                for (p1, p2) in itertools.product(points1[point1], points2[point2], repeat=1):
                    iters += 1
                    if equal_length:
                        h = hamming(p1[1], p2[1])
                    else:
                        h = 1000000
                    if h <= k:
                        reduced += 1
                        close_pairs.append((p1[0], p2[0], h))
                    elif h <= 3 * k:
                        d = lev_dis(p1[1], p2[1])
                        comp += 1
                        if d <= k:
                            close_pairs.append((p1[0], p2[0], d))
            else:
                throw += len(points1[point1]) * len(points2[point2])
    if debug:
        print("points iters  = " + str(iters))
        print("points comps  = " + str(comp))
        print("throw         = " + str(throw))
        print("len           = " + str(len(close_pairs)))
    return close_pairs


def equal_kmer_method(sequences, kmer_dict, k):
    close_pairs = []
    comps = 0
    pairs_to_compare = set()
    for name, seq in sequences.items():
        length = len(kmer_dict[name])
        possible_seq = {}
        tuples_to_sort = [(i, len(kmer_dict[kmer_dict[name][i]])) for i in range(length)]
        tuples_to_sort = sorted(tuples_to_sort, key=lambda t: t[1])
        for i in range(length):
            if i <= k:
                for possible_name in kmer_dict[kmer_dict[name][tuples_to_sort[i][0]]]:
                    if possible_name not in possible_seq:
                        possible_seq[possible_name] = 0
                    possible_seq[possible_name] += 1
            else:
                possible_seq = {key: v for key, v in possible_seq.items() if
                                i - v <= k or key in kmer_dict[kmer_dict[name][tuples_to_sort[i][0]]]}
                possible_seq = {key: v + 1 if key in kmer_dict[kmer_dict[name][tuples_to_sort[i][0]]] else v for key, v
                                in possible_seq.items()}
        for s in possible_seq:
            if name != s and (s, name) not in pairs_to_compare:
                pairs_to_compare.add((name, s))
    for n1, n2 in pairs_to_compare:
        comps += 1
        d = lev_dis(sequences[n1], sequences[n2])
        if d <= k:
            close_pairs.append((n1, n2))
    if debug:
        print("comps = " + str(comps))
        print("len = " + str(len(close_pairs)))
    return close_pairs


def equal_kmer_method_two_samples(sequences1, sequences2, kmer_dict1, kmer_dict2, points1, points2, k):
    close_pairs = []
    comps = 0
    pairs_to_compare = set()

    k_mer_coincidences = 0
    removed = 0
    for hashes_set in kmer_dict1["whole_sample"]:
        for hash_value in hashes_set:
            if hash_value in kmer_dict2:
                k_mer_coincidences += 1
                break
    if k_mer_coincidences < len(kmer_dict1["whole_sample"]) - k:
        print("k_mer_coincidences = " + str(k_mer_coincidences) + " required = " + str(
            len(kmer_dict1["whole_sample"]) - k))
        return []
    if False:
        return points_method_two_samples(sequences1, sequences2, points1, points2, k)
    sequences1 = filter_unlikly_sequences(k, kmer_dict1, kmer_dict2, sequences1)
    sequences2 = filter_unlikly_sequences(k, kmer_dict2, kmer_dict1, sequences2)
    if len(sequences1) * len(sequences2) == 0:
        return []
    if len(sequences1) < len(sequences2):
        # swap to optimize filtering
        sequences1, sequences2 = sequences2, sequences1
    kmer_dict1 = build_kmers_dict(sequences1, 14)
    kmer_dict2 = build_kmers_dict(sequences2, 14)
    for name, seq in sequences1.items():
        length = len(kmer_dict1[name])
        absence_count = 0
        # for i in range(length):
        #     if kmer_dict1[name][i] not in kmer_dict2:
        #         absence_count += 1
        # if absence_count > k:
        #     # print("absence_count "+str(absence_count)+" for "+name)
        #     removed += 1
        #     continue
        possible_seq = {}
        tuples_to_sort = [(i, len(kmer_dict2[kmer_dict1[name][i]])) for i in range(length) if
                          kmer_dict1[name][i] in kmer_dict2]
        tuples_to_sort = sorted(tuples_to_sort, key=lambda t: t[1])
        tuples_length = len(tuples_to_sort)
        sample_length = len(kmer_dict1["whole_sample"])
        for i in range(tuples_length):
            if i <= tuples_length - (sample_length - k):
                for possible_name in kmer_dict2[kmer_dict1[name][tuples_to_sort[i][0]]]:
                    if possible_name not in possible_seq:
                        possible_seq[possible_name] = 0
                    possible_seq[possible_name] += 1
            else:
                #filter
                qwe = 1
                possible_seq = {key: v for key, v in possible_seq.items() if
                                sample_length - k <= v + tuples_length - i or key in kmer_dict2[kmer_dict1[name][tuples_to_sort[i][0]]]}
                #update
                qwe = 1
                possible_seq = {key: v + 1 if key in kmer_dict2[kmer_dict1[name][tuples_to_sort[i][0]]] else v for
                                key, v in possible_seq.items()}
        for s in possible_seq:
            if name != s and (s, name) not in pairs_to_compare and possible_seq[s] >= length - k:
                pairs_to_compare.add((name, s))
    for n1, n2 in pairs_to_compare:
        comps += 1
        d = lev_dis(sequences1[n1], sequences2[n2])
        if d <= k:
            close_pairs.append((n1, n2))
    if debug:
        print("removed = " + str(removed) + " out of " + str(len(sequences1)))
        print("comps = " + str(comps))
        print("len = " + str(len(close_pairs)))
    return close_pairs


def filter_unlikly_sequences(k, kmer_dict1, kmer_dict2, sequences1):
    to_remove = set()
    for name, seq in sequences1.items():
        length = len(kmer_dict1[name])
        absence_count = 0
        for i in range(length):
            if kmer_dict1[name][i] not in kmer_dict2:
                absence_count += 1
        if absence_count > k:
            to_remove.add(name)
    return {key: v for key, v in sequences1.items() if key not in to_remove}
