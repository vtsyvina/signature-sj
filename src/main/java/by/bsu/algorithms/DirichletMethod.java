package by.bsu.algorithms;

import by.bsu.model.KMerDict;
import by.bsu.model.Sample;
import by.bsu.model.Pair;
import by.bsu.util.HammingDistance;
import by.bsu.util.KMerDictBuilder;
import by.bsu.util.LevenshteinDistance;
import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.cursors.IntCursor;
import com.carrotsearch.hppc.cursors.IntIntCursor;
import com.carrotsearch.hppc.cursors.LongCursor;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Algorithm use Dirichlet method to filter pairs on sequences from sample/samples
 */
public class DirichletMethod {

    public static boolean DEBUG = false;

    public static List<Pair> run(Sample sample1, Sample sample2, KMerDict dict1, KMerDict dict2, int k) {
        List<Pair> result = new ArrayList<>();
        int comps = 0;
        Set<Pair> pairsToCompare = new HashSet<>();
        long start = System.currentTimeMillis();
        int kMerCoincidences = 0;
        for (LongSet positionHashes : dict1.wholeSampleFixedPositionHashesList) {
            for (LongCursor hash : positionHashes) {
                if (dict2.allHashesSet.contains(hash.value)) {
                    kMerCoincidences++;
                    break;
                }
            }
        }

        if (kMerCoincidences < dict1.fixedkMersCount - k) {
            return result;
        }
        int oldLength = sample1.sequences.size() * sample2.sequences.size();
        sample1 = filterUnlikelySequences(k, dict1, dict2, sample1);
        sample2 = filterUnlikelySequences(k, dict2, dict1, sample2);

        if (sample1.sequences.size() * sample2.sequences.size() == 0) {
            return result;
        }
        if (sample1.sequences.size() < sample2.sequences.size()) {
            Sample tmp = sample1;
            sample1 = sample2;
            sample2 = tmp;
        }
        if (sample1.sequences.size() * sample2.sequences.size() < oldLength) {
            dict1 = KMerDictBuilder.getDict(sample1, dict1.l);
            dict2 = KMerDictBuilder.getDict(sample2, dict2.l);
        }
        if (DEBUG) {
            System.out.println("rebuild " + (System.currentTimeMillis() - start));
        }
        boolean sameLength = sample1.sequences.values().iterator().next().length() ==
                sample2.sequences.values().iterator().next().length();
        long construct = 0;
        long co_short = 0;
        for (Map.Entry<Integer, String> seqEntity : sample1.sequences.entrySet()) {
            IntIntMap possibleSequences = new IntIntHashMap(dict2.sequencesNumber);
            int seq = seqEntity.getKey();
            List<Pair> tuplesToSort = new ArrayList<>();
            for (int i = 0; i < dict1.fixedkMersCount; i++) {
                long hash = dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
                if (dict2.allHashesSet.contains(hash)) {
                    tuplesToSort.add(new Pair(i, dict2.hashToSequencesMap.get(hash).size()));
                }
            }
            tuplesToSort.sort(Comparator.comparing(o -> o.r));
            for (int i = 0; i < tuplesToSort.size(); i++) {
                if (i <= tuplesToSort.size() - (dict1.fixedkMersCount - k)) {
                    for (IntCursor possibleSeq :
                            dict2.hashToSequencesMap
                                    .get(dict1.sequenceFixedPositionHashesList.get(seq)[tuplesToSort.get(i).l])
                            ) {
                        possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                    }
                } else {
                    long s = System.currentTimeMillis();
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    for (IntIntCursor entry : possibleSequences) {
                        long hash = dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())[tuplesToSort.get(i).l];
                        boolean isInSecondDIct = dict2.hashToSequencesMap.get(hash).contains(entry.key);
                        if (isInSecondDIct ||
                                dict1.fixedkMersCount - k <= entry.value + tuplesToSort.size() - i) {
                            int add = isInSecondDIct ? 1 : 0;
                            tmp.put(entry.key, entry.value + add);
                        }
                    }
                    possibleSequences = tmp;
                    construct += System.currentTimeMillis() - s;
                }
            }
            for (IntIntCursor s : possibleSequences) {
                if (!seqEntity.getKey().equals(s.key)
                        && !pairsToCompare.contains(new Pair(seqEntity.getKey(), s.key))
                        && s.value >= dict1.fixedkMersCount - k) {
                    pairsToCompare.add(new Pair(seqEntity.getKey(), s.key));
                }
            }
        }
        if (DEBUG) {
            System.out.println("i < k " + co_short);
            System.out.println("i > k " + construct);
            System.out.println("constructing " + (System.currentTimeMillis() - start));
        }
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        int reduce = 0;
        for (Pair pair : pairsToCompare) {
            comps++;
            if (sameLength && hammingDistance.apply(sample1.sequences.get(pair.l), sample2.sequences.get(pair.r)) <= k) {
                result.add(pair);
                reduce++;
                continue;
            }
            int d = distance.apply(sample1.sequences.get(pair.l), sample2.sequences.get(pair.r));
            if (d != -1) {
                result.add(pair);
            }
        }
        if (DEBUG) {
            System.out.println("comparasions " + (System.currentTimeMillis() - start));
            System.out.println("comps = " + comps);
            System.out.println("reduce = " + reduce);
            System.out.println("length = " + result.size());
        }
        if (!result.isEmpty()) {
            System.out.printf("Found %s %s%n", sample1.name, sample2.name);
        }
        return result;
    }

    public static Set<Pair> run(Sample sample, KMerDict dict, int k) {
        Set<Pair> result = new HashSet<>();
        Set<Pair> pairsToCompare = new HashSet<>();
        AtomicInteger comps = new AtomicInteger(0);
        int iter = 0;
        for (Map.Entry<Integer, String> seqEntity : sample.sequences.entrySet()) {
            iter++;
            if (iter % 1000 == 0) {
                System.out.println(iter);
            }
            IntIntMap possibleSequences = new IntIntHashMap(dict.sequencesNumber);
            int seq = seqEntity.getKey();
            List<Pair> tuplesToSort = new ArrayList<>();
            for (int i = 0; i < dict.fixedkMersCount; i++) {
                long hash = dict.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
                if (dict.allHashesSet.contains(hash)) {
                    tuplesToSort.add(new Pair(i, dict.hashToSequencesMap.get(hash).size()));
                }
            }
            tuplesToSort.sort(Comparator.comparing(o -> o.r));
            for (int i = 0; i < tuplesToSort.size(); i++) {
                if (i <= tuplesToSort.size() - (dict.fixedkMersCount - k)) {
                    for (IntCursor possibleSeq :
                            dict.hashToSequencesMap
                                    .get(dict.sequenceFixedPositionHashesList.get(seq)[tuplesToSort.get(i).l])
                            ) {
                        //avoid equal pairs
                        if (possibleSeq.value == seq) {
                            continue;
                        }
                        possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    for (IntIntCursor entry : possibleSequences) {
                        long hash = dict.sequenceFixedPositionHashesList.get(seqEntity.getKey())[tuplesToSort.get(i).l];
                        boolean isInDIct = dict.hashToSequencesMap.get(hash).contains(entry.key);
                        if (isInDIct ||
                                dict.fixedkMersCount - k <= entry.value + tuplesToSort.size() - i) {
                            int add = isInDIct ? 1 : 0;
                            tmp.put(entry.key, entry.value + add);
                        }
                    }
                    possibleSequences = tmp;

                }
            }
            for (IntIntCursor s : possibleSequences) {
                if (!seqEntity.getKey().equals(s.key)
                        && !pairsToCompare.contains(new Pair(seqEntity.getKey(), s.key))
                        && s.value >= dict.fixedkMersCount - k) {
                    pairsToCompare.add(new Pair(seqEntity.getKey(), s.key));
                }
            }
        }
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        AtomicInteger reduce = new AtomicInteger(0);
        iter = 0;
        System.out.println("pairs " + pairsToCompare.size());

        pairsToCompare.stream().forEach(pair -> {

            if (hammingDistance.apply(sample.sequences.get(pair.l), sample.sequences.get(pair.r)) <= k) {
                result.add(pair);
            } else {
                comps.incrementAndGet();
                if (comps.get() % 100000 == 0) {
                    System.out.println(comps.get());
                }
                int d = distance.apply(sample.sequences.get(pair.l), sample.sequences.get(pair.r));
                if (d != -1) {
                    result.add(pair);
                }
            }
        });
        if (DEBUG) {
            System.out.println("reduce = " + reduce);
        }
        if (!result.isEmpty()) {
            System.out.printf("Found %s%n", sample.name);
            System.out.println("comps = " + comps);
            System.out.println("length = " + result.size());
        }
        return result;
    }

    public static Set<Pair> runParallel(Sample sample, KMerDict dict, int k) {
        Set<Pair> result = ConcurrentHashMap.newKeySet();
        Set<Pair> pairsToCompare = ConcurrentHashMap.newKeySet();
        AtomicInteger comps = new AtomicInteger(0);

        sample.sequences.entrySet().parallelStream().forEach(seqEntity -> {
            IntIntMap possibleSequences = new IntIntHashMap(dict.sequencesNumber);
            int seq = seqEntity.getKey();
            List<Pair> tuplesToSort = new ArrayList<>();
            for (int i = 0; i < dict.fixedkMersCount; i++) {
                long hash = dict.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
                if (dict.allHashesSet.contains(hash)) {
                    tuplesToSort.add(new Pair(i, dict.hashToSequencesMap.get(hash).size()));
                }
            }
            tuplesToSort.sort(Comparator.comparing(o -> o.r));
            for (int i = 0; i < tuplesToSort.size(); i++) {
                if (i <= tuplesToSort.size() - (dict.fixedkMersCount - k)) {
                    for (IntCursor possibleSeq :
                            dict.hashToSequencesMap
                                    .get(dict.sequenceFixedPositionHashesList.get(seq)[tuplesToSort.get(i).l])
                            ) {
                        //avoid equal pairs
                        if (possibleSeq.value == seq) {
                            continue;
                        }
                        possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    for (IntIntCursor entry : possibleSequences) {
                        long hash = dict.sequenceFixedPositionHashesList.get(seqEntity.getKey())[tuplesToSort.get(i).l];
                        boolean isInDIct = dict.hashToSequencesMap.get(hash).contains(entry.key);
                        if (isInDIct ||
                                dict.fixedkMersCount - k <= entry.value + tuplesToSort.size() - i) {
                            int add = isInDIct ? 1 : 0;
                            tmp.put(entry.key, entry.value + add);
                        }
                    }
                    possibleSequences = tmp;

                }
            }
            for (IntIntCursor s : possibleSequences) {
                if (!seqEntity.getKey().equals(s.key)
                        && !pairsToCompare.contains(new Pair(seqEntity.getKey(), s.key))
                        && s.value >= dict.fixedkMersCount - k) {
                    pairsToCompare.add(new Pair(seqEntity.getKey(), s.key));
                }
            }
        });
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        AtomicInteger reduce = new AtomicInteger(0);
        pairsToCompare.parallelStream().forEach(pair -> {
            comps.incrementAndGet();
            if (hammingDistance.apply(sample.sequences.get(pair.l), sample.sequences.get(pair.r)) <= k) {
                result.add(pair);
                reduce.incrementAndGet();

            } else {
                int d = distance.apply(sample.sequences.get(pair.l), sample.sequences.get(pair.r));
                if (d != -1) {
                    result.add(pair);
                }
            }
        });

        if (DEBUG) {
            System.out.println("reduce = " + reduce);
        }
        if (!result.isEmpty()) {
            System.out.printf("Found %s%n", sample.name);
            System.out.println("comps = " + comps);
            System.out.println("length = " + result.size());
        }
        return result;
    }

    private static Sample filterUnlikelySequences(int k, KMerDict dict1, KMerDict dict2, Sample sample) {
        Sample result = new Sample();
        result.name = sample.name;
        result.sequences = new HashMap<>();
        for (Map.Entry<Integer, String> seqEntity : sample.sequences.entrySet()) {
            int absenceCount = 0;
            for (long hash : dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())) {
                if (!dict2.hashToSequencesMap.containsKey(hash)) {
                    absenceCount++;
                }
            }
            if (absenceCount <= k) {
                result.sequences.put(seqEntity.getKey(), seqEntity.getValue());
            }
        }
        return result;
    }
}
