package by.bsu.algorithms;

import by.bsu.model.IntIntPair;
import by.bsu.model.KMerDict;
import by.bsu.model.Sample;
import by.bsu.start.Start;
import by.bsu.distance.HammingDistance;
import by.bsu.distance.LevenshteinDistance;

import com.carrotsearch.hppc.IntIntHashMap;
import com.carrotsearch.hppc.IntIntMap;
import com.carrotsearch.hppc.IntScatterSet;
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
import java.util.concurrent.atomic.AtomicInteger;

import static by.bsu.util.Utils.expandNumbers;
import static by.bsu.util.Utils.numbers;

/**
 * Algorithm use Signature method to filter pairs on sequences from sample/samples
 */
public class SignatureMethod {

    public static boolean DEBUG = false;
    private static volatile long tasksIteration = 0;
    public static AtomicInteger coincidenceFilter = new AtomicInteger();
    public static AtomicInteger executionCount = new AtomicInteger();

    public static long run(Sample sample1, Sample sample2, KMerDict dict1, KMerDict dict2, int k) throws IOException {
        int comps = 0;
        int kMerCoincidences = calculateCoincidences(dict1, dict2);
        executionCount.incrementAndGet();
        if (kMerCoincidences < dict1.fixedkMersCount - k) {
            coincidenceFilter.incrementAndGet();
            return 0;
        }
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        Path path = Start.getOutputFilename(sample1, sample2, "signature");
        StringBuilder str = new StringBuilder();
        expandNumbers(sample1.sequences.length);
        expandNumbers(sample2.sequences.length);
        int reduce = 0;
        int iteration = 0;
        int length = 0;
        for (int seq = 0; seq < sample1.sequences.length; seq++) {
            IntIntMap possibleSequences = new IntIntHashMap();
            List<IntIntPair> chunks = getSortedChunksTwoSamples(dict1, dict2, seq);
            IntSet toCompare = new IntScatterSet();
            for (int i = 0; i < chunks.size(); i++) {
                long chunkHash = dict1.sequenceFixedPositionHashesList[seq][chunks.get(i).l];
                if (i <= chunks.size() - (dict1.fixedkMersCount - k)) {
                    for (IntCursor possibleSeq : dict2.hashToSequencesMap.get(chunkHash)) {
                        possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    IntSet sequencesWithChunks = dict2.hashToSequencesMap.get(chunkHash);
                    for (IntIntCursor entry : possibleSequences) {
                        boolean isInSecondDict = sequencesWithChunks.contains(entry.key);
                        if (isInSecondDict ||
                                dict1.fixedkMersCount - k <= entry.value + chunks.size() - i) {
                            int add = isInSecondDict ? 1 : 0;
                            if (entry.value + add >= dict1.fixedkMersCount - k) {
                                toCompare.add(entry.key);
                            } else {
                                tmp.put(entry.key, entry.value + add);
                            }
                        }
                    }
                    possibleSequences = tmp;
                }
            }
            for (IntIntCursor s : possibleSequences) {
                if (seq != s.key
                        && s.value >= dict1.fixedkMersCount - k) {
                    comps++;
                    if (hammingDistance.apply(sample1.forHamming[seq], sample2.forHamming[s.key]) <= k) {
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                        reduce++;
                        length++;
                        continue;
                    }
                    int d = distance.apply(sample1.sequences[seq], sample2.sequences[s.key]);
                    if (d != -1) {
                        length++;
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                    }
                }
            }
            iteration++;
            if (iteration % 400 == 0) {
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
            }
        }
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        if (DEBUG && comps > 0) {
            System.out.printf("%s %s%n", sample1.name, sample2.name);
            System.out.println("comps = " + comps);
            System.out.println("reduce = " + reduce);
            System.out.println("levenshtein = " + (comps - reduce));
            System.out.println("length = " + length);
        }
        if (length > 0) {
            System.out.printf("Found %s %s. Length = %d\n", sample1.name, sample2.name, length);
        } else {
            Files.delete(path);
        }
        return length;
    }

    public static long run(Sample sample, KMerDict dict, int k) throws IOException {
        System.out.println("Start Signature method for " + sample.name + " k=" + k + " l=" + dict.l);
        expandNumbers(sample.sequences.length);
        long[] iter = {0, 0, 0, 0};
        int[] distances = new int[264];
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        StringBuilder str = new StringBuilder();
        Path path = Start.getOutputFilename(sample, "signature");
        long length = 0;
        iter[0] = 0;
        for (int seq = 0; seq < sample.sequences.length; seq++) {
            iter[0]++;
            //write to file each 400 iterations
            if (iter[0] % 400 == 0) {
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r" + iter[0]);
            }
            IntIntMap possibleSequences = new IntIntHashMap();
            List<IntIntPair> sortedChunks = getSortedChunksOneSample(dict, seq);
            IntSet toCompare = new IntScatterSet();
            for (int i = 0; i < sortedChunks.size(); i++) {
                if (i <= sortedChunks.size() - (dict.fixedkMersCount - k)) {
                    fillPossiblePairs(dict, possibleSequences, seq, sortedChunks, i);
                } else {
                    possibleSequences = filterPossibleSequences(dict, possibleSequences, seq, k, sortedChunks, i, toCompare);
                }
            }
            String s1 = sample.sequences[seq];
            String h1 = sample.forHamming[seq];
            for (IntCursor s : toCompare) {
                iter[1]++;
                if (hammingDistance.apply(h1, sample.forHamming[s.value]) <= k) {
                    length++;
                    iter[2]++;
                    str.append(numbers.get(seq)).append(" ").append(numbers.get(s.value)).append("\n");
                } else {
                    if (distance.apply(s1, sample.sequences[s.value]) != -1) {
                        length++;
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.value)).append("\n");
                    } else if (DEBUG) {
                        //int z = unlim.apply(sample.sequences.get(seq), sample.sequences.get(s.key));
                        //distances[z]++;
                    }
                }
            }
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
            System.out.println("levenshtein = " + (iter[1] - iter[2] - iter[3]));
            System.out.println("c = " + iter[3]);
        }
        return length;
    }

    public static Long runParallel(Sample sample, KMerDict dict, int k) throws IOException {
        System.out.println("Start Signature method parallel for " + sample.name + " k= " + k + " l= " + dict.l);
        expandNumbers(sample.sequences.length);
        Path path = Start.getOutputFilename(sample, "signature");
        //divide sequences into parts for executor service
        int cores = Runtime.getRuntime().availableProcessors();
        ExecutorService service = Executors.newFixedThreadPool(cores);
        List<Map<Integer, String>> parts = new ArrayList<>();
        tasksIteration = 0;
        int partsCount = Math.min(cores, sample.sequences.length / 400 + 1);
        for (int i = 0; i < partsCount; i++) {
            parts.add(new HashMap<>());
        }
        for (int i = 0; i < sample.sequences.length; i++) {
            parts.get(i % partsCount).put(i, sample.sequences[i]);
        }
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
            System.out.println("levenshtein = " + (results[1] - results[2]));
            System.out.println("length = " + results[3]);
        }
        return results[0];
    }

    /**
     * Fills possibleSequences with appearing sequences for given position(chunks.get(iter))
     */
    private static void fillPossiblePairs(KMerDict dict, IntIntMap possibleSequences, int seq, List<IntIntPair> chunks, int iter) {
        for (IntCursor possibleSeq :
                dict.hashToSequencesMap
                        .get(dict.sequenceFixedPositionHashesList[seq][chunks.get(iter).l])
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
     * <p>
     * add sequence to toCompare set if it already reached necessary amount of hits so we don't need it in possibleSequences anymore
     */
    private static IntIntMap filterPossibleSequences(KMerDict dict, IntIntMap possibleSequences, int seq, int k, List<IntIntPair> chunks, int iter, IntSet toCompare) {
        IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
        long hash = dict.sequenceFixedPositionHashesList[seq][chunks.get(iter).l];
        IntSet sequencesWithHashSet = dict.hashToSequencesMap.get(hash);
        for (IntIntCursor entry : possibleSequences) {
            boolean isInDict = sequencesWithHashSet.contains(entry.key);
            // put if sequence hash l-mer for current fixed position or if it already has enough equal l-mers
            if (isInDict ||
                    dict.fixedkMersCount - k <= entry.value + chunks.size() - iter) {
                int add = isInDict ? 1 : 0;
                if (entry.value + add >= dict.fixedkMersCount - k) {
                    toCompare.add(entry.key);
                } else {
                    tmp.put(entry.key, entry.value + add);
                }
            }
        }
        return tmp;
    }

    /**
     * Returns list of sorted chunks for current sequence entity by number of sequences that contain
     * l-mers from sequence entity fixed positions
     */
    private static List<IntIntPair> getSortedChunksOneSample(KMerDict dict, int seq) {
        List<IntIntPair> result = new ArrayList<>();
        //for each fixed position
        for (int i = 0; i < dict.fixedkMersCount; i++) {
            long hash = dict.sequenceFixedPositionHashesList[seq][i];
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
     * see getSortedChunksOneSample
     */
    private static List<IntIntPair> getSortedChunksTwoSamples(KMerDict dict1, KMerDict dict2, int seq) {
        List<IntIntPair> result = new ArrayList<>();
        for (int i = 0; i < dict1.fixedkMersCount; i++) {
            long hash = dict1.sequenceFixedPositionHashesList[seq][i];
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
                List<IntIntPair> sortedChunks = getSortedChunksOneSample(dict, seq);
                IntSet toCompare = new IntScatterSet();
                for (int i = 0; i < sortedChunks.size(); i++) {
                    if (i <= sortedChunks.size() - (dict.fixedkMersCount - k)) {
                        fillPossiblePairs(dict, possibleSequences, seq, sortedChunks, i);
                    } else {
                        possibleSequences = filterPossibleSequences(dict, possibleSequences, seq, k, sortedChunks, i, toCompare);
                    }
                }
                String s1 = sample.sequences[seq];
                String h1 = sample.forHamming[seq];
                for (IntCursor s : toCompare) {
                    iters[1]++;
                    if (hammingDistance.apply(h1, sample.forHamming[s.value]) <= k) {
                        iters[3]++;
                        iters[2]++;
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.value)).append("\n");
                    } else {
                        if (distance.apply(s1, sample.sequences[s.value]) != -1) {
                            iters[3]++;
                            str.append(numbers.get(seq)).append(" ").append(numbers.get(s.value)).append("\n");
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
