package by.bsu.algorithms;

import static by.bsu.util.AlgorithmUtils.calculateCoincidences;
import static by.bsu.util.Utils.expandNumbers;
import static by.bsu.util.Utils.numbers;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import com.carrotsearch.hppc.IntIntHashMap;
import com.carrotsearch.hppc.IntIntMap;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.cursors.IntIntCursor;

import by.bsu.distance.HammingDistance;
import by.bsu.model.IntIntPair;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;
import by.bsu.start.Start;

/**
 * Created by c5239200 on 6/26/17.
 * Algorithm based on the same idea that in signature method,
 * but find only sequences that are closed in terms of Hamming Distance
 */
public class SignatureHammingMethod {

    public static boolean DEBUG = false;
    private static volatile long tasksIteration = 0;
    public static AtomicInteger coincidenceFilter = new AtomicInteger();
    public static AtomicInteger executionCount = new AtomicInteger();
    
    //TODO check this method. Just copy-past for now, and docs for whole class
    public static long run(Sample sample1, Sample sample2, KMerDictChunks dict1, KMerDictChunks dict2, int k) throws IOException {
        int comps = 0;
        int kMerCoincidences = calculateCoincidences(dict1, dict2);
        executionCount.incrementAndGet();
        if (kMerCoincidences < dict1.fixedkMersCount - k) {
            coincidenceFilter.incrementAndGet();
            return 0;
        }
        HammingDistance hammingDistance = new HammingDistance();
        Path path = Start.getOutputFilename(sample1, sample2, "signature");
        StringBuilder str = new StringBuilder();
        expandNumbers(sample1.sequences.length);
        expandNumbers(sample2.sequences.length);
        int iteration = 0;
        int length = 0;
        int[] hits = new int[sample2.sequences.length];
        for (int seq = 0; seq < sample1.sequences.length; seq++) {
            IntIntMap possibleSequences = new IntIntHashMap();
            List<IntIntPair> chunks = getSortedChunksTwoSamples(dict1, dict2, seq);
            for (int i = 0; i < chunks.size(); i++) {
                long chunkHash = dict1.sequenceFixedPositionHashesList[seq][chunks.get(i).l];
                if (i <= chunks.size() - (dict1.fixedkMersCount - k)) {
                    for (int possibleSeq : dict2.chunksHashToSequencesMapArray[chunks.get(i).l].get(chunkHash)) {
                        hits[possibleSeq]++;
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    int[] sequencesWithChunks = dict2.chunksHashToSequencesMap[chunks.get(i).l].get(chunkHash);
                    for (IntIntCursor entry : possibleSequences) {
                        boolean isInSecondDict = sequencesWithChunks[entry.key] == 1;
                        if (isInSecondDict ||
                                dict1.fixedkMersCount - k <= entry.value + chunks.size() - i) {
                            int add = isInSecondDict ? 1 : 0;
                            tmp.put(entry.key, entry.value + add);
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
                        length++;
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
            System.out.println("length = " + length);
        }
        if (length > 0) {
            System.out.printf("Found %s %s. Length = %d\n", sample1.name, sample2.name, length);
        } else {
            Files.delete(path);
        }
        return length;
    }

    public static long run(Sample sample, KMerDictChunks dict, int k) throws IOException {
        System.out.println("Start Signature Hamming method for " + sample.name + " k=" + k + " l=" + dict.l);
        expandNumbers(sample.sequences.length);
        long[] iter = {0, 0};
        int[] distances = new int[dict.sequencesLength];
        HammingDistance hammingDistance = new HammingDistance();
        StringBuilder str = new StringBuilder();
        Path path = Start.getOutputFilename(sample, "signature-hamming");
        long length = 0;
        iter[0] = 0;
        long hammingTime = 0;
        long fillTime = 0;
        long filterTime = 0;
        long writeTime = 0;
        for (int seq = 0; seq < sample.sequences.length ; seq++) {
            iter[0]++;
            //write to file each 400 iterations
            if (iter[0] % 400 == 0) {
                long start = System.nanoTime();
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r" + iter[0]);
                writeTime += System.nanoTime() - start;
            }
            int[] possibleSequences = new int[1];
            int[] hits = new int[sample.sequences.length];
            List<IntIntPair> sortedChunks = getSortedChunksOneSample(dict, seq);
            List<Integer> toCompare = new ArrayList<>();
            for (int i = 0; i < sortedChunks.size(); i++) {
                if (i <= sortedChunks.size() - (dict.fixedkMersCount - k)) {
                    long start = System.nanoTime();
                    fillPossiblePairs(dict, hits, seq, sortedChunks, i);
                    fillTime += System.nanoTime() - start;
                } else {
                    long start = System.nanoTime();
                    if (i-1 == sortedChunks.size() - (dict.fixedkMersCount - k)){
                        possibleSequences = fillPossibleSequencesFromArray(hits);
                    }
                    possibleSequences = filterPossibleSequences(dict, possibleSequences, hits, seq, k, sortedChunks, i, toCompare);
                    filterTime += System.nanoTime() - start;
                }
            }
            //String h1 = sample.forHamming.get(seq);
            char[] h1 = sample.sequencesChars[seq];
            long[] c1 = dict.sequenceFixedPositionHashesList[seq];
            long start = System.currentTimeMillis();
            for (Integer s : toCompare) {
                iter[1]++;
                int apply = hammingDistance.apply(h1, sample.sequencesChars[s],
                        c1, dict.sequenceFixedPositionHashesList[s],
                        dict.l, k);
                if (apply != -1) {
                    length++;
                    str.append(numbers.get(seq)).append(" ").append(numbers.get(s)).append("\n");
                }
            }
            hammingTime += System.currentTimeMillis() - start;
        }
        //write the rest of computed pairs
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println();
        if (length != 0) {
            System.out.printf("Found %s%n", sample.name);
            System.out.println("comps = " + iter[1]);
            System.out.println("length = " + length);
            System.out.println("hamming time = " + hammingTime);
            System.out.println("fill time = " + fillTime/1000000);
            System.out.println("filter time = " + filterTime/1000000);
            System.out.println("write time = " + writeTime/1000000);
            //System.out.println("distances = " + Arrays.toString(distances));
        }
        return length;
    }

    /**
     * Using array instead of collection provide performance improvement due to low cost of iteration and random access by index
     */
    private static int[] fillPossibleSequencesFromArray(int[] fill) {
        int[] tmp = new int[fill.length];
        int last = 0;
        for (int j = 0; j < fill.length; j++) {
            if (fill[j] != 0){
                tmp[last++] = j;
            }
        }
        int[] result = new int[last];
        System.arraycopy(tmp, 0, result, 0, last);
        return result;
    }

    public static Long runParallel(Sample sample, KMerDictChunks dict, int k) throws IOException {
        System.out.println("Start Signature Hamming method parallel for " + sample.name + " k= " + k + " l= " + dict.l);
        expandNumbers(sample.sequences.length);
        Path path = Start.getOutputFilename(sample, "signature-hamming");
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
        parts.forEach(part -> tasks.add(new SignatureHammingMethod.ParallelTask(sample, part, dict, k, path)));
        long[] results = {0, 0, 0, 0, 0};
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
            System.out.println("length = " + results[3]);
        }
        return results[0];
    }

    /**
     * Increments possibleSequences array with appearing sequences for given position(chunks.get(iter))
     */
    private static void fillPossiblePairs(KMerDictChunks dict, int[] hits, int seq, List<IntIntPair> chunks, int iter) {
        for (int possibleSeq :
                dict.chunksHashToSequencesMapArray[chunks.get(iter).l]
                        .get(dict.sequenceFixedPositionHashesList[seq][chunks.get(iter).l])
                ) {
            //avoid equal pairs
            if (possibleSeq <= seq) {
                continue;
            }
            hits[possibleSeq]++;
        }
    }

    /**
     * Filters possibleSequences by removing all sequences that doesn't contain necessary amount of equal l-mers
     * with current sequence (seq)
     * <p>
     * add sequence to toCompare set if it already reached necessary amount of hits so we don't need it in possibleSequences anymore
     */
    private static int[] filterPossibleSequences(KMerDictChunks dict, int[] possibleSequences, int[] hits, int seq, int k, List<IntIntPair> chunks, int iter, List<Integer> toCompare) {
        int[] tmp = new int[possibleSequences.length];
        long hash = dict.sequenceFixedPositionHashesList[seq][chunks.get(iter).l];
        int[] sequencesWithHashSet = dict.chunksHashToSequencesMap[chunks.get(iter).l].get(hash);
        int last = 0;
        for (int candidate : possibleSequences) {
            boolean isInDict = sequencesWithHashSet[candidate] == 1;
            // put if sequence hash l-mer for current fixed position or if it already has enough equal l-mers
            if (isInDict ||
                    dict.fixedkMersCount - k <= hits[candidate] + chunks.size() - iter) {
                if (isInDict){
                    hits[candidate]++;
                }
                if (hits[candidate] >= dict.fixedkMersCount - k) {
                    toCompare.add(candidate);
                } else {
                    tmp[last++] = candidate;
                }
            }
        }
        int[] result = new int[last];
        System.arraycopy(tmp, 0, result, 0, last);
        return result;
    }

    /**
     * Returns list of sorted chunks for current sequence entity by number of sequences that contain
     * l-mers from sequence entity fixed positions
     */
    private static List<IntIntPair> getSortedChunksOneSample(KMerDictChunks dict, int seq) {
        List<IntIntPair> result = new ArrayList<>();
        //for each fixed position
        for (int i = 0; i < dict.fixedkMersCount; i++) {
            long hash = dict.sequenceFixedPositionHashesList[seq][i];
            if (dict.allHashesSet.contains(hash)) {
                //add tuple -> (position, amount of sequences)
                result.add(new IntIntPair(i, dict.chunksHashToSequencesMapArray[i].get(hash).length));
            }
        }
        //sort by amount
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * see getSortedChunksOneSample
     */
    private static List<IntIntPair> getSortedChunksTwoSamples(KMerDictChunks dict1, KMerDictChunks dict2, int seq) {
        List<IntIntPair> result = new ArrayList<>();
        for (int i = 0; i < dict1.fixedkMersCount; i++) {
            long hash = dict1.sequenceFixedPositionHashesList[seq][i];
            if (dict2.allHashesSet.contains(hash)) {
                result.add(new IntIntPair(i, dict2.chunksHashToSequencesMapArray[i].get(hash).length));
            }
        }
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * Class runs almost the same code as sequential run, but on different sequences set
     */
    private static class ParallelTask implements Callable<long[]> {
        private Sample sample;
        private Map<Integer, String> sequences;
        private KMerDictChunks dict;
        private int k;
        private final Path path;

        ParallelTask(Sample sample, Map<Integer, String> sequences, KMerDictChunks dict, int k, Path path) {
            this.sample = sample;
            this.sequences = sequences;
            this.dict = dict;
            this.k = k;
            this.path = path;
        }

        @Override
        public long[] call() throws Exception {
            HammingDistance hammingDistance = new HammingDistance();
            StringBuilder str = new StringBuilder();
            /*
              0 -> iteration
              1 -> comparisons
              2 -> hamming
              3 -> total length
             */
            long[] iters = {0, 0, 0, 0, 0};
            //only for additional QGram filter experiment(doesn't work good)
            //int q = 8;
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
                int[] possibleSequences = new int[1];
                int seq = seqEntity.getKey();
                List<IntIntPair> sortedChunks = getSortedChunksOneSample(dict, seq);
                List<Integer> toCompare = new ArrayList<>();
                int[] hits = new int[sample.sequences.length];
                for (int i = 0; i < sortedChunks.size(); i++) {
                    if (i <= sortedChunks.size() - (dict.fixedkMersCount - k)) {
                        fillPossiblePairs(dict, hits, seq, sortedChunks, i);
                    } else {
                        if (i-1 == sortedChunks.size() - (dict.fixedkMersCount - k)){
                            possibleSequences = fillPossibleSequencesFromArray(hits);
                        }
                        possibleSequences = filterPossibleSequences(dict, possibleSequences, hits, seq, k, sortedChunks, i, toCompare);
                    }
                }
                //String h1 = sample.forHamming.get(seq);
                char[] h1 = sample.sequencesChars[seq];
                long[] c1 = dict.sequenceFixedPositionHashesList[seq];
                for (Integer s : toCompare) {
                    iters[1]++;
                    int apply = hammingDistance.apply(h1, sample.sequencesChars[s],
                            c1, dict.sequenceFixedPositionHashesList[s],
                            dict.l, k);
                    if (apply != -1) {
                        iters[3]++;
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s)).append("\n");
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