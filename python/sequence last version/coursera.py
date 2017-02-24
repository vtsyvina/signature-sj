import random


# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i + k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates


# Input:  A list Items
# Output: A list containing all objects from Items without duplicates
def remove_duplicates(Items):
    ItemsNoDuplicates = set(Items)  # output variable
    return ItemsNoDuplicates


# Input:  A string Text and an integer k
# Output: CountDict(Text, k)
# HINT:   This code should be identical to when you last implemented CountDict
def CountDict(Text, k):
    Count = {}
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i + k]
        Count[i] = PatternCount(Pattern, Text)
    return Count


# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
# HINT:   This code should be identical to when you last implemented PatternCount
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count = count + 1
    return count


def PatternMatching(Pattern, Genome):
    positions = []  # output variable
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i + (n // 2)])
    return array


def HammingDistance(p, q):
    result = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            result += 1
    return result


# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Pattern, Text, d):
    positions = []  # initializing list of positions
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i + len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions


def MinimumSkew(Genome):
    positions = []  # output variable
    skew = Skew(Genome)
    min = 100000;
    for key, value in skew.items():
        if value < min:
            min = value
            positions = [key]
        elif value == min:
            positions.append(key)
    return positions


# Input:  A String Genome
# Output: Skew(Genome)
# HINT:   This code should be taken from the last Code Challenge.
def Skew(Genome):
    skew = {}  # initializing the dictionary
    skew[0] = 0;
    for i in range(len(Genome)):
        if Genome[i] == 'G':
            skew[i + 1] = skew[i] + 1
        elif Genome[i] == 'C':
            skew[i + 1] = skew[i] - 1
        else:
            skew[i + 1] = skew[i]
    return skew


def Count(Motifs):
    count = {}  # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] += 1
    return profile


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    score = 0
    for i in range(t):
        for j in range(k):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score


# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] /= t
    return profile


# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}  # output variable
    profile = CountWithPseudocounts(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] /= t + 4
    return profile


# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    # insert your code here
    count = Count(Motifs)
    consensus = ""
    k = len(Motifs[0])
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    t = len(Text)
    pr = 1
    for i in range(t):
        pr *= Profile[Text[i]][i]
    return pr


# Input:  String Text, an integer k, and profile matrix Profile
# Output: ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    maxPr = -1
    result = ""
    for i in range(len(Text) - k + 1):
        if Pr(Text[i:(i + k)], Profile) > maxPr:
            maxPr = Pr(Text[i:(i + k)], Profile)
            result = Text[i:(i + k)]
    return result


# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []  # output variable
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    result = []
    n = len(Dna)
    for i in range(n):
        result.append(ProfileMostProbablePattern(Dna[i], len(Profile['A']), Profile))
    return result


def RandomMotifs(Dna, k, t):
    result = []
    l = len(Dna[0])
    for i in range(t):
        left = random.randint(0, l - k)
        result.append(Dna[i][left:left + k])
    return result


# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    sum = 0
    result = {}
    for kmer in Probabilities:
        sum += Probabilities[kmer]
    for kmer in Probabilities:
        result[kmer] = Probabilities[kmer] / sum
    return result


# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = ''  # output variable
    norm = Normalize(Probabilities)
    sum = 0
    r = random.uniform(0, 1)
    for k in norm:
        if sum <= r < sum + norm[k]:
            kmer = k
            break
        sum += norm[k]
    return kmer


# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    BestMotifs = []  # output variable
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(1, N):
        i = random.randint(1, t - 1)
        profile = Profile(Motifs[:i] + Motifs[i + 1:])
        Motifs[i] = ProfileGeneratedString(Dna[i], profile, k)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
