package by.bsu.algorithms;

import static by.bsu.util.Utils.expandNumbers;
import static by.bsu.util.Utils.numbers;

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

import com.carrotsearch.hppc.IntIntHashMap;
import com.carrotsearch.hppc.IntIntMap;
import com.carrotsearch.hppc.IntScatterSet;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongSet;
import com.carrotsearch.hppc.cursors.IntCursor;
import com.carrotsearch.hppc.cursors.IntIntCursor;
import com.carrotsearch.hppc.cursors.LongCursor;

import by.bsu.model.IntIntPair;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;
import by.bsu.start.Start;
import by.bsu.distance.HammingDistance;

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
        expandNumbers(sample1.sequences.size());
        expandNumbers(sample2.sequences.size());
        int iteration = 0;
        int length = 0;
        for (Map.Entry<Integer, String> seqEntity : sample1.sequences.entrySet()) {
            IntIntMap possibleSequences = new IntIntHashMap();
            int seq = seqEntity.getKey();
            List<IntIntPair> chunks = getSortedChunksTwoSamples(dict1, dict2, seqEntity);
            for (int i = 0; i < chunks.size(); i++) {
                long chunkHash = dict1.sequenceFixedPositionHashesList.get(seq)[chunks.get(i).l];
                if (i <= chunks.size() - (dict1.fixedkMersCount - k)) {
                    for (IntCursor possibleSeq : dict2.chunksHashToSequencesMap[chunks.get(i).l].get(chunkHash)) {
                        possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    IntSet sequencesWithChunks = dict2.chunksHashToSequencesMap[chunks.get(i).l].get(chunkHash);
                    for (IntIntCursor entry : possibleSequences) {
                        boolean isInSecondDict = sequencesWithChunks.contains(entry.key);
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
                    if (hammingDistance.apply(sample1.forHamming.get(seq), sample2.forHamming.get(s.key)) <= k) {
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
        expandNumbers(sample.sequences.size());
        long[] iter = {0, 0};
        int[] distances = new int[dict.sequencesLength];
        HammingDistance hammingDistance = new HammingDistance();
        StringBuilder str = new StringBuilder();
        Path path = Start.getOutputFilename(sample, "signature-hamming");
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
            IntIntMap possibleSequences = new IntIntHashMap();
            int seq = seqEntity.getKey();
            List<IntIntPair> sortedChunks = getSortedChunksOneSample(dict, seqEntity);
            IntSet toCompare = new IntScatterSet();
            for (int i = 0; i < sortedChunks.size(); i++) {
                if (i <= sortedChunks.size() - (dict.fixedkMersCount - k)) {
                    fillPossiblePairs(dict, possibleSequences, seq, sortedChunks, i);
                } else {
                    possibleSequences = filterPossibleSequences(dict, possibleSequences, seq, k, sortedChunks, i, toCompare);
                }
            }
            //String h1 = sample.forHamming.get(seq);
            char[] h1 = sample.sequencesChars.get(seq);
            long[] c1 = dict.sequenceFixedPositionHashesList.get(seq);
            for (IntCursor s : toCompare) {
                iter[1]++;
                int apply = hammingDistance.apply(h1, sample.sequencesChars.get(s.value),
                        c1, dict.sequenceFixedPositionHashesList.get(s.value),
                        dict.l, k);
                if (apply != -1) {
                    length++;
                    str.append(numbers.get(seq)).append(" ").append(numbers.get(s.value)).append("\n");
                }
            }
        }
        //write the rest of computed pairs
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println();
        if (length != 0) {
            System.out.printf("Found %s%n", sample.name);
            System.out.println("comps = " + iter[1]);
            System.out.println("length = " + length);
            System.out.println("distances = " + Arrays.toString(distances));
        }
        return length;
    }

    public static Long runParallel(Sample sample, KMerDictChunks dict, int k) throws IOException {
        System.out.println("Start Signature Hamming method parallel for " + sample.name + " k= " + k + " l= " + dict.l);
        expandNumbers(sample.sequences.size());
        Path path = Start.getOutputFilename(sample, "signature-hamming");
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
     * Fills possibleSequences with appearing sequences for given position(chunks.get(iter))
     */
    private static void fillPossiblePairs(KMerDictChunks dict, IntIntMap possibleSequences, int seq, List<IntIntPair> chunks, int iter) {
        for (IntCursor possibleSeq :
                dict.chunksHashToSequencesMap[chunks.get(iter).l]
                        .get(dict.sequenceFixedPositionHashesList.get(seq)[chunks.get(iter).l])
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
    private static IntIntMap filterPossibleSequences(KMerDictChunks dict, IntIntMap possibleSequences, int seq, int k, List<IntIntPair> chunks, int iter, IntSet toCompare) {
        IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
        long hash = dict.sequenceFixedPositionHashesList.get(seq)[chunks.get(iter).l];
        IntSet sequencesWithHashSet = dict.chunksHashToSequencesMap[chunks.get(iter).l].get(hash);
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
    private static List<IntIntPair> getSortedChunksOneSample(KMerDictChunks dict, Map.Entry<Integer, String> seqEntity) {
        List<IntIntPair> result = new ArrayList<>();
        //for each fixed position
        for (int i = 0; i < dict.fixedkMersCount; i++) {
            long hash = dict.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
            if (dict.allHashesSet.contains(hash)) {
                //add tuple -> (position, amount of sequences)
                result.add(new IntIntPair(i, dict.chunksHashToSequencesMap[i].get(hash).size()));
            }
        }
        //sort by amount
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * see getSortedChunksOneSample
     */
    private static List<IntIntPair> getSortedChunksTwoSamples(KMerDictChunks dict1, KMerDictChunks dict2, Map.Entry<Integer, String> seqEntity) {
        List<IntIntPair> result = new ArrayList<>();
        for (int i = 0; i < dict1.fixedkMersCount; i++) {
            long hash = dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
            if (dict2.allHashesSet.contains(hash)) {
                result.add(new IntIntPair(i, dict2.chunksHashToSequencesMap[i].get(hash).size()));
            }
        }
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * Calculates how many l-mers from fixed positions from first dictionary are in second dictionary
     */
    private static int calculateCoincidences(KMerDictChunks dict1, KMerDictChunks dict2) {
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
                IntIntMap possibleSequences = new IntIntHashMap();
                int seq = seqEntity.getKey();
                List<IntIntPair> sortedChunks = getSortedChunksOneSample(dict, seqEntity);
                IntSet toCompare = new IntScatterSet();
                for (int i = 0; i < sortedChunks.size(); i++) {
                    if (i <= sortedChunks.size() - (dict.fixedkMersCount - k)) {
                        fillPossiblePairs(dict, possibleSequences, seq, sortedChunks, i);
                    } else {
                        possibleSequences = filterPossibleSequences(dict, possibleSequences, seq, k, sortedChunks, i, toCompare);
                    }
                }
                //String h1 = sample.forHamming.get(seq);
                char[] h1 = sample.sequencesChars.get(seq);
                long[] c1 = dict.sequenceFixedPositionHashesList.get(seq);
                for (IntCursor s : toCompare) {
                    iters[1]++;
                    int apply = hammingDistance.apply(h1, sample.sequencesChars.get(s.value),
                            c1, dict.sequenceFixedPositionHashesList.get(s.value),
                            dict.l, k);
                    if (apply != -1) {
                        iters[3]++;
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.value)).append("\n");
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
