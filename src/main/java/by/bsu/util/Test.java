package by.bsu.util;

import by.bsu.model.Sample;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import static java.util.stream.Collectors.toMap;

/**
 * Created by c5239200 on 2/3/17.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@State(Scope.Thread)
@Warmup(iterations = 10, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 20, time = 500, timeUnit = TimeUnit.MILLISECONDS)
public class Test {
    LevenshteinDistance defaultLev = LevenshteinDistance.getDefaultInstance();
    LevenshteinDistance limitedLev = new LevenshteinDistance(30);
    HammingDistance hammingDistance = new HammingDistance();
    String str1 = "CACCGACTGCGGCGCTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGAACATCGCGCTTATCA";
    String str2 = "CACCGACTGCGGCACTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGGACATCGCGCTTATCA";
    AtomicInteger i = new AtomicInteger(0);
    Sample query = null;
    long[] h1;
    long[] h2;
    {
        Map<Integer, String > seq;
        Map<Integer, String > tmp;
        try {
            seq = FasReader.readList(Paths.get("test_data/db8/32000.fas"));
            tmp = FasReader.readList(Paths.get("test_data/db8/32000 (2).fas"));
            int[] size = new int[1];
            size[0] = seq.size();
            tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
            seq.putAll(tmp);
            tmp = FasReader.readList(Paths.get("test_data/db8/32000 (3).fas"));
            size[0] = seq.size();
            tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
            seq.putAll(tmp);
            tmp = FasReader.readList(Paths.get("test_data/db8/32000 (4).fas"));
            size[0] = seq.size();
            tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
            seq.putAll(tmp);
            query = new Sample("db8", seq);
        } catch (IOException e) {
            e.printStackTrace();
        }
        h1 = Utils.getSequenceHashesArray(query.sequences.get(1), 11);
        h2 = Utils.getSequenceHashesArray(query.sequences.get(2), 11);
    }

    //@Benchmark
    @Fork(1)
    public void testLevenshtein() {
        defaultLev.apply(str1, str2);
    }

    //@Benchmark
    @Fork(1)
    public void testLevenshteinLimited() {
        limitedLev.apply(str1, str2);
    }

    //@Benchmark
    @Fork(1)
    public void testHammingLimited() {
        hammingDistance.apply(str1, str2);
    }


    //@Benchmark
    @Fork(1)
    public void testAtomic() {
        i.incrementAndGet();
    }

    @Benchmark
    @Fork(1)
    public void testConsensus() {
        Utils.distancesMap(query.consensus, query.sequences, 20);
    }

    @Benchmark
    @Fork(1)
    public void testConsensus2() {
        Utils.distancesMap(query.consensus, query.sequences, 20);
    }
}
