import itertools

from Levenshtein import hamming
from Levenshtein import distance as lev_dis
from functions import lower_bound_estimate_list, sequence_coordinates, lower_bound_estimate, build_segments_dict

debug = True

def points_method_additional_coords_union(sequences, k):
    points = {}
    for name, sequence in sequences.items():
        full_coords = get_full_coords(sequence)
        if full_coords not in points:
            points[full_coords] = []

        points[full_coords].append(
            (name, sequence, full_coords))
    close_pairs = []
    iters = 0
    comp = 0

    to_delete = []
    to_append = {}
    for point in points:
        if len(points[point]) == 1:
            temp = (point[0], (0, 0, 0, 0), (0, 0, 0, 0))
            if temp not in to_append:
                to_append[temp] = []
            to_append[temp].append(points[point][0])
            to_delete.append(point)
    for point in to_delete:
        del points[point]
    for point in to_append:
        points[point] = to_append[point]
    count = {}
    n = 0
    for point in points:
        if len(points[point]) not in count:
            count[len(points[point])] = 0
        count[len(points[point])] += 1
        n += len(points[point])
    to_process = set(points)
    for point in points:
        for point2 in to_process:
            iters += 1
            if point == point2:
                for (p1, p2) in itertools.combinations(points[point], r=2):
                    iters += 1
                    if hamming(p1[1], p2[1]) <= 2 * k:
                        d = lev_dis(p1[1], p2[1])
                        comp += 1
                        if d <= k:
                            close_pairs.append((p1[0], p2[0], d))
            elif lower_bound_estimate_list(point[0], point2[0]) <= k \
                    and (point[1] == (0, 0, 0, 0) or point[2] == (0, 0, 0, 0)
                        or point2[1] == (0, 0, 0, 0) or point2[2] == (0, 0, 0, 0)
                         or (lower_bound_estimate_list(point[1], point2[1]) <= k
                             and lower_bound_estimate_list(point[2], point2[2]) <= k)):
                for (p1, p2) in itertools.product(points[point], points[point2], repeat=1):
                    iters += 1
                    if True or hamming(p1[1], p2[1]) <= 2 * k:
                        d = lev_dis(p1[1], p2[1])
                        comp += 1
                        if d <= k:
                            close_pairs.append((p1[0], p2[0], d))
                        pass
        to_process.remove(point)
    if debug:
        print("points iters  = " + str(iters))
        print("points comps  = " + str(comp))
        print(len(close_pairs))
    return close_pairs

def points_method_additional_coords(sequences, k):
    points = {}
    for name, sequence in sequences.items():
        full_coords = get_full_coords(sequence)
        if full_coords not in points:
            points[full_coords] = []

        points[full_coords].append(
            (name, sequence, full_coords))
    close_pairs = []
    iters = 0
    comp = 0
    count = {}
    n = 0
    for point in points:
        if len(points[point]) not in count:
            count[len(points[point])] = 0
        count[len(points[point])] += 1
        n += len(points[point])
    if debug:
        print("percent of all " + str(count[1] / n))
    to_process = set(points)
    for point in points:
        for point2 in to_process:
            iters += 1
            if point == point2:
                for (p1, p2) in itertools.combinations(points[point], r=2):
                    iters += 1
                    if hamming(p1[1], p2[1]) <= 2 * k:
                        d = lev_dis(p1[1], p2[1])
                        comp += 1
                        if d <= k:
                            close_pairs.append((p1[0], p2[0], d))
            elif lower_bound_estimate_list(point[0], point2[0]) <= k \
                    and lower_bound_estimate_list(point[1], point2[1]) <= k \
                    and lower_bound_estimate_list(point[2], point2[2]) <= k:
                for (p1, p2) in itertools.product(points[point], points[point2], repeat=1):
                    iters += 1
                    if hamming(p1[1], p2[1]) <= 2 * k:
                        d = lev_dis(p1[1], p2[1])
                        comp += 1
                        if d <= k:
                            close_pairs.append((p1[0], p2[0], d))
        to_process.remove(point)
    if debug:
        print("points iters  = " + str(iters))
        print("points comps  = " + str(comp))
        print(len(close_pairs))
    return close_pairs

def neighbourhood_method(sequences, k):
    neighbourhoods = {}
    for name, sequence in sequences.items():
        coors = sequence_coordinates(sequence)
        vertex = k + 1
        for (x1, x2, x3, x4) in itertools.product([0, vertex], repeat=4):

            y1 = (coors['A'] // vertex) * vertex + x1
            y2 = (coors['C'] // vertex) * vertex + x2
            y3 = (coors['G'] // vertex) * vertex + x3
            y4 = (coors['T'] // vertex) * vertex + x4
            if not (y1, y2, y3, y4) in neighbourhoods:
                neighbourhoods[(y1, y2, y3, y4)] = []
            if abs(y1 - coors['A']) <= k \
                    and abs(y2 - coors['C']) <= k \
                    and abs(y3 - coors['G']) <= k \
                    and abs(y4 - coors['T']) <= k:
                neighbourhoods[(y1, y2, y3, y4)].append((name, coors))
    if debug:
        print("neighbourhoods len " + str(len(neighbourhoods)))
    pairs_to_check = set()
    iters_to_pair_check = 0
    for coors, secuences_list in neighbourhoods.items():
        for ((name1, coord1), (name2, coord2)) in itertools.combinations(secuences_list, 2):
            iters_to_pair_check += 1
            if lower_bound_estimate(coord1, coord2) <= k \
                    and (name1, name2) not in pairs_to_check \
                    and (name2, name1) not in pairs_to_check \
                    and hamming(sequences[name1], sequences[name2]) <= 2 * k:
                pairs_to_check.add((name1, name2))
    if debug:
        print("iters_to_pair_check " + str(iters_to_pair_check))
        print("pairs_to_check len " + str(len(pairs_to_check)))
    close_pairs = []
    comp = 0
    for (name1, name2) in pairs_to_check:

        d = lev_dis(sequences[name1], sequences[name2])
        comp += 1
        if d <= k:
            close_pairs.append(((name1, sequences[name1]), (name2, sequences[name2]), d))
    if debug:
        print("ne comp = " + str(comp))
        print(len(close_pairs))
    return close_pairs

def get_full_coords(sequence):
    coords = sequence_coordinates(sequence)
    coords = (coords['A'], coords['C'], coords['G'], coords['T'])
    coords_left = sequence_coordinates(sequence[len(sequence) // 2:])
    coords_left = (coords_left['A'], coords_left['C'], coords_left['G'], coords_left['T'])
    coords_right = sequence_coordinates(sequence[:len(sequence) // 2])
    coords_right = (coords_right['A'], coords_right['C'], coords_right['G'], coords_right['T'])
    full_coords = (coords, coords_left, coords_right)
    return full_coords

def common_k_mer_method_two_samples(sequences1, sequences2, k):
    piece_length = 20
    pieces_count = len(list(sequences1.values())[0]) // piece_length
    pairs_to_compare = {}
    result = []
    for pieces_combination in itertools.combinations( range(pieces_count), r= pieces_count-k):
        d1 = build_segments_dict(sequences1, pieces_combination[0], piece_length)
        d2 = build_segments_dict(sequences2, pieces_combination[0], piece_length)
        result.extend(recursive_k_mer_move(d1,d2, 1, pieces_combination, piece_length))
    # print(result)
    for (d1,d2) in result:
        for (name1, seq1) in d1.items():
            for (name2, seq2) in d2.items():
                if (name1,name2) not in pairs_to_compare:
                    pairs_to_compare[(name1,name2)] = (seq1, seq2)
    # print(pairs_to_compare)

def recursive_k_mer_move(d1, d2, step, combination, piece_length):
    if step == len(combination):
        (min_dict, max_dict) = (d1, d2) if len(d1) < len(d2) else (d2, d1)
        result = []
        for segment in min_dict:
            if segment in max_dict:
                result.append((min_dict[segment], max_dict[segment]))
        return result
    (min_dict, max_dict) = (d1,d2) if len(d1) < len(d2) else (d2,d1)
    result = []
    for segment in min_dict:
        if segment in max_dict:
            next_step_d1 = build_segments_dict(d1[segment], combination[step], piece_length)
            next_step_d2 = build_segments_dict(d2[segment], combination[step], piece_length)
            s = []
            if (len(next_step_d1) > 0 and len(next_step_d2) > 0):
                s = recursive_k_mer_move(next_step_d1, next_step_d2, step+1, combination, piece_length)
            if len(s) > 0:
                result.extend(s)
    return result