package by.bsu.algorithms;

import by.bsu.distance.HammingDistance;
import by.bsu.model.IntIntPair;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;
import by.bsu.start.Start;
import by.bsu.util.AlgorithmUtils;
import com.carrotsearch.hppc.IntIntHashMap;
import com.carrotsearch.hppc.IntIntMap;
import com.carrotsearch.hppc.cursors.IntIntCursor;

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

import static by.bsu.util.AlgorithmUtils.calculateCoincidences;
import static by.bsu.util.Utils.expandNumbers;
import static by.bsu.util.Utils.numbers;

/**
 * Created by c5239200 on 6/26/17.
 * Algorithm based on the same idea that in signature method,
 * but find only sequences that are closed in terms of Hamming Distance
 */
public class SignatureHammingMethod {

    private Path output;
    private Sample sample1;
    private Sample sample2;

    public static boolean DEBUG = false;
    private static volatile long tasksIteration = 0;
    public static AtomicInteger coincidenceFilter = new AtomicInteger();
    public static AtomicInteger executionCount = new AtomicInteger();

    public long run(Sample sample1, Sample sample2, KMerDictChunks dict1, KMerDictChunks dict2, int k) throws IOException {
        this.sample1 = sample1;
        this.sample2 = sample2;
        int comps = 0;
        int kMerCoincidences = calculateCoincidences(dict1, dict2);
        executionCount.incrementAndGet();
        if (kMerCoincidences < dict1.chunksCount - k) {
            coincidenceFilter.incrementAndGet();
            return 0;
        }
        HammingDistance hammingDistance = new HammingDistance();
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
                long chunkHash = dict1.sequenceChunksHashesList[seq][chunks.get(i).l];
                if (i <= chunks.size() - (dict1.chunksCount - k)) {
                    for (int possibleSeq : dict2.chunksHashToSequencesMapArray[chunks.get(i).l].get(chunkHash)) {
                        hits[possibleSeq]++;
                    }
                } else {
                    IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                    byte[] sequencesWithChunks = dict2.chunksHashToSequencesMap[chunks.get(i).l].get(chunkHash);
                    for (IntIntCursor entry : possibleSequences) {
                        boolean isInSecondDict = sequencesWithChunks[entry.key] == 1;
                        if (isInSecondDict ||
                                dict1.chunksCount - k <= entry.value + chunks.size() - i) {
                            int add = isInSecondDict ? 1 : 0;
                            tmp.put(entry.key, entry.value + add);
                        }
                    }
                    possibleSequences = tmp;
                }
            }
            for (IntIntCursor s : possibleSequences) {
                if (seq != s.key
                        && s.value >= dict1.chunksCount - k) {
                    comps++;
                    if (hammingDistance.apply(sample1.forHamming[seq], sample2.forHamming[s.key]) <= k) {
                        str.append(numbers.get(seq)).append(" ").append(numbers.get(s.key)).append("\n");
                        length++;
                    }
                }
            }
            iteration++;
            if (iteration % 400 == 0) {
                Files.write(getPath(), str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
            }
        }
        Files.write(getPath(), str.toString().getBytes(), StandardOpenOption.APPEND);
        if (DEBUG && comps > 0) {
            System.out.printf("%s %s%n", sample1.name, sample2.name);
            System.out.println("comps = " + comps);
            System.out.println("related pairs found = " + length);
        }
        if (length > 0) {
            System.out.printf("Found %s %s. Length = %d\n", sample1.name, sample2.name, length);
        } else {
            Files.delete(getPath());
        }
        return length;
    }

    public long run(Sample sample, KMerDictChunks dict, int k) throws IOException {
        String chunks = dict.l == 0 ? " entropy-based chunks size" : " l=" + dict.l;
        System.out.println("Start Signature Hamming method for " + sample.name + " k=" + k + chunks);
        System.out.println("Input size = " + sample.sequences.length);
        expandNumbers(sample.sequences.length);
        long[] iter = {0, 0};
        HammingDistance hammingDistance = new HammingDistance();
        StringBuilder str = new StringBuilder();
        Path path = Start.getOutputFilename(sample, "signature-hamming");
        long length = 0;
        iter[0] = 0;
        long hammingTime = 0;
        long fillTime = 0;
        long filterTime = 0;
        long writeTime = 0;
        for (int seq = 0; seq < sample.sequences.length; seq++) {
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
                if (i <= sortedChunks.size() - (dict.chunksCount - k)) {
                    long start = System.nanoTime();
                    fillPossiblePairs(dict, hits, seq, sortedChunks, i);
                    fillTime += System.nanoTime() - start;
                } else {
                    long start = System.nanoTime();
                    if (i - 1 == sortedChunks.size() - (dict.chunksCount - k)) {
                        possibleSequences = AlgorithmUtils.fillPossibleSequencesFromArray(hits);
                    }
                    possibleSequences = filterPossibleSequences(dict, possibleSequences, hits, seq, k, sortedChunks, i, toCompare);
                    filterTime += System.nanoTime() - start;
                }
            }
            long start = System.currentTimeMillis();
            for (Integer s : toCompare) {
                iter[1]++;
                int apply = hammingDistance.apply(sample.forHamming[seq], sample.forHamming[s]);
                if (apply <= k) {
                    length++;
                    str.append(numbers.get(seq)).append(" ").append(numbers.get(s)).append("\n");
                }
            }
            hammingTime += System.currentTimeMillis() - start;
        }
        //write the rest of computed pairs
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println();
        System.out.println("comps = " + iter[1]);
        System.out.println("related pairs found = " + length);
//            System.out.println("hamming time = " + hammingTime);
//            System.out.println("fill time = " + fillTime / 1000000);
//            System.out.println("filter time = " + filterTime / 1000000);
//            System.out.println("write time = " + writeTime / 1000000);
        System.out.println("Output is available at " + path.toAbsolutePath().toString());
        return length;
    }

    public Long runParallel(Sample sample, KMerDictChunks dict, int k) throws IOException {
        String chunks = dict.l == 0 ? " entropy-based segments size" : " l=" + dict.l;
        System.out.println("Start Signature Hamming method for " + sample.name + " k=" + k + chunks);
        System.out.println("Input size = " + sample.sequences.length);
        expandNumbers(sample.sequences.length);
        Path path = Start.getOutputFilename(sample, "signature-hamming");
        int cores = Runtime.getRuntime().availableProcessors();

        String threads = Start.settings.get("-threads");
        if (threads != null) {
            cores = Integer.valueOf(threads);
        }
        System.out.println("Running threads = " + cores);
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

        System.out.println("comparisons = " + results[1]);
        System.out.println("related pairs found = " + results[3]);
        System.out.println("Output is available at " + path.toAbsolutePath().toString());
        return results[0];
    }

    /**
     * Increments possibleSequences array with appearing sequences for given position(chunks.get(iter))
     */
    private void fillPossiblePairs(KMerDictChunks dict, int[] hits, int seq, List<IntIntPair> chunks, int iter) {
        for (int possibleSeq :
                dict.chunksHashToSequencesMapArray[chunks.get(iter).l]
                        .get(dict.sequenceChunksHashesList[seq][chunks.get(iter).l])
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
    private int[] filterPossibleSequences(KMerDictChunks dict, int[] possibleSequences, int[] hits, int seq, int k, List<IntIntPair> chunks, int iter, List<Integer> toCompare) {
        int[] tmp = new int[possibleSequences.length];
        long hash = dict.sequenceChunksHashesList[seq][chunks.get(iter).l];
        byte[] sequencesWithHashSet = dict.chunksHashToSequencesMap[chunks.get(iter).l].get(hash);
        int last = 0;
        for (int candidate : possibleSequences) {
            boolean isInDict = sequencesWithHashSet[candidate] == 1;
            // put if sequence hash l-mer for current fixed position or if it already has enough equal l-mers
            if (isInDict ||
                    dict.chunksCount - k <= hits[candidate] + chunks.size() - iter) {
                if (isInDict) {
                    hits[candidate]++;
                }
                if (hits[candidate] >= dict.chunksCount - k) {
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
    private List<IntIntPair> getSortedChunksOneSample(KMerDictChunks dict, int seq) {
        List<IntIntPair> result = new ArrayList<>();
        //for each fixed position
        for (int i = 0; i < dict.chunksCount; i++) {
            long hash = dict.sequenceChunksHashesList[seq][i];
            //add tuple -> (position, amount of sequences)
            result.add(new IntIntPair(i, dict.chunksHashToSequencesMapArray[i].get(hash).length));
        }
        //sort by amount
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * see getSortedChunksOneSample
     */
    private List<IntIntPair> getSortedChunksTwoSamples(KMerDictChunks dict1, KMerDictChunks dict2, int seq) {
        List<IntIntPair> result = new ArrayList<>();
        for (int i = 0; i < dict1.chunksCount; i++) {
            long hash = dict1.sequenceChunksHashesList[seq][i];
            if (dict2.chunksHashToSequencesMapArray[i].containsKey(hash)) {
                result.add(new IntIntPair(i, dict2.chunksHashToSequencesMapArray[i].get(hash).length));
            }
        }
        result.sort(Comparator.comparing(o -> o.r));
        return result;
    }

    /**
     * Class runs almost the same code as sequential run, but on different sequences set
     */
    private class ParallelTask implements Callable<long[]> {
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
                    if (i <= sortedChunks.size() - (dict.chunksCount - k)) {
                        fillPossiblePairs(dict, hits, seq, sortedChunks, i);
                    } else {
                        if (i - 1 == sortedChunks.size() - (dict.chunksCount - k)) {
                            possibleSequences = AlgorithmUtils.fillPossibleSequencesFromArray(hits);
                        }
                        possibleSequences = filterPossibleSequences(dict, possibleSequences, hits, seq, k, sortedChunks, i, toCompare);
                    }
                }
                //String h1 = sample.forHamming.get(seq);
                for (Integer s : toCompare) {
                    iters[1]++;
                    int apply = hammingDistance.apply(sample.forHamming[seq], sample.forHamming[s]);
                    if (apply <= k) {
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

    private Path getPath() throws IOException {
        if (output == null) {
            output = Start.getOutputFilename(sample1, sample2, "signature-hamming");
        }
        return output;
    }
}
