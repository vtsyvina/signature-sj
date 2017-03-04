package by.bsu.util;

import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by c5239200 on 2/3/17.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@State(Scope.Thread)
@Warmup(iterations = 10, time = 200, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 20, time = 200, timeUnit = TimeUnit.MILLISECONDS)
public class Test {
    LevenshteinDistance defaultLev = LevenshteinDistance.getDefaultInstance();
    LevenshteinDistance limitedLev = new LevenshteinDistance(0);
    HammingDistance hammingDistance = new HammingDistance();
    String str1 = "CACCGACTGCGGCGCTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGAACATCGCGCTTATCA";
    String str2 = "CACCGACTGCGGCACTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGGACATCGCGCTTATCA";
    AtomicInteger i = new AtomicInteger(0);
    @Benchmark
    @Fork(1)
    public void testLevenshtein() {
        defaultLev.apply(str1, str2);
    }

    @Benchmark
    @Fork(1)
    public void testLevenshteinLimited() {
        limitedLev.apply(str1, str2);
    }

    @Benchmark
    @Fork(1)
    public void testHammingLimited() {
        hammingDistance.apply(str1, str2);
    }


    @Benchmark
    @Fork(1)
    public void testAtomic() {
        i.incrementAndGet();
    }
}
