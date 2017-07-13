package by.bsu.util;

import by.bsu.algorithms.QGramSimilarity;
import by.bsu.distance.HammingDistance;
import by.bsu.distance.LevenshteinDistance;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;
import info.debatty.java.stringsimilarity.QGram;
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

import com.sun.scenario.effect.impl.sw.sse.SSEBlend_SRC_OUTPeer;

/**
 * Created by c5239200 on 2/3/17.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@State(Scope.Thread)
@Warmup(iterations = 5, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Measurement(iterations = 10, time = 500, timeUnit = TimeUnit.MILLISECONDS)
public class Test {
    LevenshteinDistance defaultLev = LevenshteinDistance.getDefaultInstance();
    LevenshteinDistance limitedLev = new LevenshteinDistance(11);
    HammingDistance hammingDistance = new HammingDistance();

    QGram qGram7 = new QGram(7);
    QGramSimilarity qGramSimilarity = new QGramSimilarity();
    String str1 = "CACCGACTGCGGCGCTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGAACATCGCGCTTATCA";
    String str2 = "CACCGACTGCGGCACTGGTTATGGCACAAGTGCTCCGGATCCCGGAAGCTATCGTGGATATGGTAGCTGGAGCCCACTGGGGAGTCCTAGCGGGGCTAGCTTACTATTCCATGGTTGGCAACTGGGCGAAGGTGCTAGTCGTGCTGCTCCTGTTCGCGGGGGTTGATGCTGATACCAAGACCATCGGCGGTAAGGCTACGCAGCAAACCGCGCGCCTCACCAGCTTCTTTAGCCCGGGTCCCCAGCAGGACATCGCGCTTATCA";
    Map<String, Integer> p1;
    Map<String, Integer> p2;
    AtomicInteger i = new AtomicInteger(0);
    Sample query = null;
    KMerDictChunks dict;
    long[] h1;
    long[] h2;
    char[] s1 = new char[264];
    char[] s2 = new char[264];
    int l = 11;
    int k = 4;
    int left = 1;
    int right = 279;
    {
        p1 = qGram7.getProfile(str1);
        p2 = qGram7.getProfile(str2);
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
            dict = KMerDictChunksBuilder.getDict(query, l);
            str1 = query.sequences.get(left);
            str2 = query.sequences.get(right);
            str1.getChars(0, str1.length(), s1, 0);
            str2.getChars(0, str2.length(), s2, 0);
        } catch (IOException e) {
            e.printStackTrace();
        }
        h1 = dict.sequenceFixedPositionHashesList.get(left);
        h2 = dict.sequenceFixedPositionHashesList.get(right);
    }

    //@Benchmark
    @Fork(1)
    public void testLevenshtein() {
        defaultLev.apply(str1, str2);
    }
    

    //@Benchmark
    @Fork(1)
    public void testQram11() {
        qGram7.distance(p1, p2);
    }

    //@Benchmark
    @Fork(1)
    public void testQramSimilarity() {
        qGramSimilarity.similarity(p1, p2);
    }

    //@Benchmark
    @Fork(1)
    public void testLevenshteinLimited() {
        limitedLev.apply(str1, str2);
    }

    @Benchmark
    @Fork(1)
    public void testHamming() {
        hammingDistance.apply(str1, str2);
    }

    @Benchmark
    @Fork(1)
    public void testHammingWithDict() {
        long[] chunk1 = h1;
        long[] chunk2 = h2;
        str1 = query.sequences.get(left);
        str2 = query.sequences.get(right);
        int c = 0;
        it:for (int i1 = 0; i1 < chunk1.length; i1++) {
            if (chunk1[i1] != chunk2[i1]){
                c++;
                for (int j = 0; j < l; j++) {
                    if (str1.charAt(i1*l +j) != str2.charAt(i1*l + j)){
                        c++;
                        if (c > k){
                            break it;
                        }
                    }
                }
            }
        }
    }

    @Benchmark
    @Fork(1)
    public void testHammingWithDictCharArray() {
        hammingDistance.apply(s1, s2, h1, h2, l, k);
    }


    //@Benchmark
    @Fork(1)
    public void testAtomic() {
        i.incrementAndGet();
    }

    //@Benchmark
    @Fork(1)
    public void testConsensus() {
        Utils.distancesMap(query.consensus, query.sequences, 20);
    }

    //@Benchmark
    @Fork(1)
    public void testConsensus2() {
        Utils.distancesMap(query.consensus, query.sequences, 20);
    }
}
