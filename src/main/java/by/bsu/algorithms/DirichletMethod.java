package by.bsu.algorithms;

import by.bsu.model.IntIntPair;
import by.bsu.model.KMerDict;
import by.bsu.model.Sample;
import by.bsu.start.Start;
import by.bsu.util.HammingDistance;
import by.bsu.util.KMerDictBuilder;
import by.bsu.util.LevenshteinDistance;
import com.carrotsearch.hppc.IntIntHashMap;
import com.carrotsearch.hppc.IntIntMap;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongSet;
import com.carrotsearch.hppc.cursors.IntCursor;
import com.carrotsearch.hppc.cursors.IntIntCursor;
import com.carrotsearch.hppc.cursors.LongCursor;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import static by.bsu.util.Utils.expandNumbers;
import static by.bsu.util.Utils.numbers;

/**
 * Algorithm use Dirichlet method to filter pairs on sequences from sample/samples
 */
public class DirichletMethod {

    public static boolean DEBUG = true;
    private static volatile long tasksIteration = 0;

    public static Set<IntIntPair> run(Sample sample1, Sample sample2, KMerDict dict1, KMerDict dict2, int k) {
        Set<IntIntPair> result = new HashSet<>();
        int comps = 0;
        int kMerCoincidences = calculateCoincidences(dict1, dict2);
        if (kMerCoincidences < dict1.fixedkMersCount - k) {
            return result;
        }
        int oldLength = sample1.sequences.size() * sample2.sequences.size();
        sample1 = filterUnlikelySequences(k, dict1, dict2, sample1);
        sample2 = filterUnlikelySequences(k, dict2, dict1, sample2);

        //if we filter one of samples fully
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
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        int reduce = 0;
        for (Map.Entry<Integer, String> seqEntity : sample1.sequences.entrySet()) {
            IntIntMap possibleSequences = new IntIntHashMap(dict2.sequencesNumber);
            int seq = seqEntity.getKey();
            List<IntIntPair> tuples = getSortedTuplesTwoSamples(dict1, dict2, seqEntity);
            for (int i = 0; i < tuples.size(); i++) {
                if (i <= tuples.size() - (dict1.fixedkMersCount - k)) {
                    for (IntCursor possibleSeq :
                            dict2.hashToSequencesMap
                                    .get(dict1.sequenceFixedPositionHashesList.get(seq)[tuples.get(i).l])
                            ) {
                        possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    long hash = dict1.sequenceFixedPositionHashesList.get(seq)[tuples.get(i).l];
                    for (IntIntCursor entry : possibleSequences) {
                        boolean isInSecondDict = dict2.hashToSequencesMap.get(hash).contains(entry.key);
                        if (isInSecondDict ||
                                dict1.fixedkMersCount - k <= entry.value + tuples.size() - i) {
                            int add = isInSecondDict ? 1 : 0;
                            tmp.put(entry.key, entry.value + add);
                        }
                    }
                    possibleSequences = tmp;
                }
            }
            for (IntIntCursor s : possibleSequences) {
                if (seq !=s.key
                        && s.value >= dict1.fixedkMersCount - k) {
                    comps++;
                    if ( hammingDistance.apply(sample1.forHamming.get(seq), sample2.forHamming.get(s.key)) <= k) {
                        result.add(new IntIntPair(seq, s.key));
                        reduce++;
                        continue;
                    }
                    int d = distance.apply(sample1.sequences.get(seq), sample2.sequences.get(s.key));
                    if (d != -1) {
                        result.add(new IntIntPair(seq, s.key));
                    }
                }
            }
        }
        if (DEBUG && comps > 0) {
            System.out.printf("%s %s%n", sample1.name, sample2.name);
            System.out.println("comps = " + comps);
            System.out.println("reduce = " + reduce);
            System.out.println("length = " + result.size());
        }
        if (!result.isEmpty()) {
            System.out.printf("Found %s %s%n", sample1.name, sample2.name);
        }
        return result;
    }

    public static long run(Sample sample, KMerDict dict, int k) throws IOException {
        System.out.println("Start Dirihlet method for " + sample.name + " k=" + k + " l=" + dict.l);
        expandNumbers(sample.sequences.size());
        Set<Integer> processed = new HashSet<>();
        long[] iter = {0, 0, 0, 0};
        int[] distances = new int[264];
        LevenshteinDistance distance = new LevenshteinDistance(k);
        LevenshteinDistance unlim = new LevenshteinDistance(60);
        HammingDistance hammingDistance = new HammingDistance();
        StringBuilder str = new StringBuilder();
        Path path = Start.getOutputFilename(sample,"dirichlet");
        long length = 0;
        iter[0] = 0;
        for (Map.Entry<Integer, String> seqEntity : sample.sequences.entrySet()) {
            iter[0]++;
            //write to file each 400 iterations
            if (iter[0] % 400 == 0) {
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r" + iter[0]);
            }
            IntIntMap possibleSequences = new IntIntHashMap(dict.sequencesNumber);
            int seq = seqEntity.getKey();
            List<IntIntPair> sortedTuples = getSortedTuplesOneSample(dict, seqEntity);
            for (int i = 0; i < sortedTuples.size(); i++) {
                if (i <= sortedTuples.size() - (dict.fixedkMersCount - k)) {
                    fillPossiblePairs(dict, possibleSequences, seq, sortedTuples, processed, i);
                } else {
                    possibleSequences = filterPossibleSequences(dict, possibleSequences, seq, k, sortedTuples, i);
                }
            }
            String s1 = sample.sequences.get(seq);
            String h1 = sample.forHamming.get(seq);
            for (IntIntCursor s : possibleSequences) {
                if (s.value >= dict.fixedkMersCount - k) {
                    iter[1]++;
                    if (hammingDistance.apply(h1, sample.forHamming.get(s.key)) <= k) {
                        length++;
                        iter[2]++;
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                    } else {
                        if (distance.apply(s1, sample.sequences.get(s.key)) != -1) {
                            length++;
                            str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                        } else if(DEBUG) {
                            int z = unlim.apply(sample.sequences.get(seq), sample.sequences.get(s.key));
                            distances[z]++;
                        }
                    }

                }
            }
            processed.add(seq);
        }
        //write the rest of computed pairs
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println();
        if (length != 0) {
            System.out.printf("Found %s%n", sample.name);
            System.out.println("comps = " + iter[1]);
            System.out.println("reduce = " + iter[2]);
            System.out.println("length = " + length);
            System.out.println("distances = " + Arrays.toString(distances));
            System.out.println("levenshtein = " + (iter[1]-iter[2]));
        }
        return length;
    }

    public static Long runParallel(Sample sample, KMerDict dict, int k) throws IOException {
        System.out.println("Start Dirihlet method parallel for " + sample.name + " k= " + k + " l= " + dict.l);
        expandNumbers(sample.sequences.size());
        Path path = Start.getOutputFilename(sample, "dirichlet");
        //divide sequences into parts for executor service
        int cores = Runtime.getRuntime().availableProcessors();
        ExecutorService service = Executors.newFixedThreadPool(cores);
        List<Map<Integer, String>> parts = new ArrayList<>();
        tasksIteration = 0;
        int partsCount = Math.min(cores, sample.sequences.size() / 400 + 1);
        for (int i = 0; i < partsCount; i++) {
            parts.add(new HashMap<>());
        }
        final int[] i = {0};
        sample.sequences.entrySet().forEach(entity -> {
            parts.get(i[0] % partsCount).put(entity.getKey(), entity.getValue());
            i[0]++;
        });
        List<Callable<long[]>> tasks = new ArrayList<>();
        //create tasks with parts of sequences
        parts.forEach(part -> tasks.add(new ParallelTask(sample, part, dict, k, path)));
        long[] results = {0, 0, 0, 0};
        try {
            List<Future<long[]>> futures = service.invokeAll(tasks);
            service.shutdown();
            futures.forEach(future -> {
                try {
                    long[] f = future.get();
                    /*
                      0 -> iteration
                      1 -> comparisons
                      2 -> hamming
                      3 -> total length
                     */
                    results[0] += f[0];
                    results[1] += f[1];
                    results[2] += f[2];
                    results[3] += f[3];
                } catch (InterruptedException | ExecutionException e) {
                    System.err.println("Error! Parallel tasks were not successful on get");
                    e.printStackTrace();
                }
            });
        } catch (InterruptedException e) {
            System.err.println("Error! Parallel tasks were not successful on invoke");
            e.printStackTrace();
        }
        System.out.println();
        if (results[3] > 0) {
            System.out.printf("Found %s%n", sample.name);
            System.out.println("comparisons = " + results[1]);
            System.out.println("passed hamming distance = " + results[2]);
            System.out.println("levenshtein = " + (results[1]-results[2]));
            System.out.println("length = " + results[3]);
        }
        return results[0];
    }

    /**
     * Fills possibleSequences with appearing sequences for given position(tuples.get(iter))
     * only for non parallel run
     */
    private static void fillPossiblePairs(KMerDict dict, IntIntMap possibleSequences, int seq, List<IntIntPair> tuples, Set<Integer> processed, int iter) {
        for (IntCursor possibleSeq :
                dict.hashToSequencesMap
                        .get(dict.sequenceFixedPositionHashesList.get(seq)[tuples.get(iter).l])
                ) {
            //avoid equal pairs, do not add already processed sequences
            if (possibleSeq.value <= seq || processed.contains(possibleSeq.value)) {
                continue;
            }
            possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
        }
    }

    /**
     * For parallel(processed is buggy in this case)
     */
    private static void fillPossiblePairs(KMerDict dict, IntIntMap possibleSequences, int seq, List<IntIntPair> tuples, int iter) {
        for (IntCursor possibleSeq :
                dict.hashToSequencesMap
                        .get(dict.sequenceFixedPositionHashesList.get(seq)[tuples.get(iter).l])
                ) {
            //avoid equal pairs
            if (possibleSeq.value <= seq) {
                continue;
            }
            possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
        }
    }

    /**
     * Filters possibleSequences by removing all sequences that doesn't contain necessary amount of equal l-mers
     * with current sequence (seq)
     */
    private static IntIntMap filterPossibleSequences(KMerDict dict, IntIntMap possibleSequences, int seq, int k, List<IntIntPair> tuples, int iter) {
        IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
        long hash = dict.sequenceFixedPositionHashesList.get(seq)[tuples.get(iter).l];
        IntSet sequencesWithHashSet = dict.hashToSequencesMap.get(hash);
        for (IntIntCursor entry : possibleSequences) {
            boolean isInDict = sequencesWithHashSet.contains(entry.key);
            // put if sequence hash l-mer for current fixed position or if it already has enough equal l-mers
            if (isInDict ||
                    dict.fixedkMersCount - k <= entry.value + tuples.size() - iter) {
                int add = isInDict ? 1 : 0;
                tmp.put(entry.key, entry.value + add);
            }
        }
        return tmp;
    }

    /**
     * Returns list of sorted tuples for current sequence entity by number of sequences that contain
     * l-mers from sequence entity fixed positions
     */
    private static List<IntIntPair> getSortedTuplesOneSample(KMerDict dict, Map.Entry<Integer, String> seqEntity) {
        List<IntIntPair> result = new ArrayList<>();
        //for each fixed position
        for (int i = 0; i < dict.fixedkMersCount; i++) {
            long hash = dict.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
            if (dict.allHashesSet.contains(hash)) {
                //add tuple -> (position, amount of sequences)
                result.add(new IntIntPair(i, dict.hashToSequencesMap.get(hash).size()));
            }
        }
        //sort by amount
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * removes all sequences from given sample that don't have any related sequences in second sample based on
     * amount of presented l-mers for fixed positions
     */
    private static Sample filterUnlikelySequences(int k, KMerDict dict1, KMerDict dict2, Sample sample) {
        Sample result = new Sample();
        result.name = sample.name;
        result.sequences = new HashMap<>();
        result.forHamming = new HashMap<>();
        for (Map.Entry<Integer, String> seqEntity : sample.sequences.entrySet()) {
            int absenceCount = 0;
            for (long hash : dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())) {
                if (!dict2.hashToSequencesMap.containsKey(hash)) {
                    absenceCount++;
                }
            }
            if (absenceCount <= k) {
                result.sequences.put(seqEntity.getKey(), seqEntity.getValue());
                result.forHamming.put(seqEntity.getKey(), sample.forHamming.get(seqEntity.getKey()));
            }
        }
        return result;
    }

    /**
     * see getSortedTuplesOneSample
     */
    private static List<IntIntPair> getSortedTuplesTwoSamples(KMerDict dict1, KMerDict dict2, Map.Entry<Integer, String> seqEntity) {
        List<IntIntPair> result = new ArrayList<>();
        for (int i = 0; i < dict1.fixedkMersCount; i++) {
            long hash = dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
            if (dict2.allHashesSet.contains(hash)) {
                result.add(new IntIntPair(i, dict2.hashToSequencesMap.get(hash).size()));
            }
        }
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * Calculates how many l-mers from fixed positions from first dictionary are in second dictionary
     */
    private static int calculateCoincidences(KMerDict dict1, KMerDict dict2) {
        int kMerCoincidences = 0;
        for (LongSet positionHashes : dict1.wholeSampleFixedPositionHashesList) {
            for (LongCursor hash : positionHashes) {
                if (dict2.allHashesSet.contains(hash.value)) {
                    kMerCoincidences++;
                    break;
                }
            }
        }
        return kMerCoincidences;
    }



    /**
     * Class runs almost the same code as sequential run, but on different sequences set
     */
    private static class ParallelTask implements Callable<long[]> {
        private Sample sample;
        private Map<Integer, String> sequences;
        private KMerDict dict;
        private int k;
        private final Path path;

        ParallelTask(Sample sample, Map<Integer, String> sequences, KMerDict dict, int k, Path path) {
            this.sample = sample;
            this.sequences = sequences;
            this.dict = dict;
            this.k = k;
            this.path = path;
        }

        @Override
        public long[] call() throws Exception {
            LevenshteinDistance distance = new LevenshteinDistance(k);
            HammingDistance hammingDistance = new HammingDistance();
            StringBuilder str = new StringBuilder();
            /*
              0 -> iteration
              1 -> comparisons
              2 -> hamming
              3 -> total length
             */
            long[] iters = {0, 0, 0, 0};

            int fileWriteThreshold = Math.min(sequences.size() / 10, 400);
            for (Map.Entry<Integer, String> seqEntity : sequences.entrySet()) {
                iters[0]++;
                tasksIteration++;
                //write to file each fileWriteThreshold iterations
                if (iters[0] % fileWriteThreshold == 0) {
                    synchronized (path) {
                        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                    }
                    str = new StringBuilder();
                    System.out.print("\r" + tasksIteration);
                }
                IntIntMap possibleSequences = new IntIntHashMap();
                int seq = seqEntity.getKey();
                List<IntIntPair> sortedTuples = getSortedTuplesOneSample(dict, seqEntity);
                for (int i = 0; i < sortedTuples.size(); i++) {
                    if (i <= sortedTuples.size() - (dict.fixedkMersCount - k)) {
                        fillPossiblePairs(dict, possibleSequences, seq, sortedTuples, i);
                    } else {
                        possibleSequences = filterPossibleSequences(dict, possibleSequences, seq, k, sortedTuples, i);
                    }
                }
                String s1 = sample.sequences.get(seq);
                String h1 = sample.forHamming.get(seq);
                for (IntIntCursor s : possibleSequences) {
                    if (s.value >= dict.fixedkMersCount - k) {
                        iters[1]++;
                        if (hammingDistance.apply(h1, sample.forHamming.get(s.key)) <= k) {
                            iters[3]++;
                            iters[2]++;
                            str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                        } else {
                            if (distance.apply(s1, sample.sequences.get(s.key)) != -1) {
                                iters[3]++;
                                str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                            }
                        }
                    }
                }
            }
            //write the rest of computed pairs
            synchronized (path) {
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
            }
            return iters;
        }
    }

}
