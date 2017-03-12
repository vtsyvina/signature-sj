package by.bsu.algorithms;

import by.bsu.model.IntIntPair;
import by.bsu.model.Sample;
import by.bsu.util.LevenshteinDistance;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Brute force algorithm that campares all possible sequence pairs
 */
public class BruteForce {

    public static Set<IntIntPair> run(Sample sample1, Sample sample2, int k){
        Set<IntIntPair> result = new HashSet<>();
        LevenshteinDistance distance = new LevenshteinDistance(k);
        for (Map.Entry<Integer, String> entry1 : sample1.sequences.entrySet()){
            for (Map.Entry<Integer, String> entry2 : sample2.sequences.entrySet()){
                int d = distance.apply(entry1.getValue(), entry2.getValue());
                if (d != -1 && d<=k){
                    result.add(new IntIntPair(entry1.getKey(), entry2.getKey()));
                }
            }
        }
        System.out.println("length = " +result.size());
        return result;
    }

    public static Set<IntIntPair> run(Sample sample, int k){
        Set<IntIntPair> result = new HashSet<>();
        LevenshteinDistance distance = new LevenshteinDistance(k);
        for (Map.Entry<Integer, String> entry1 : sample.sequences.entrySet()){
            for (Map.Entry<Integer, String> entry2 : sample.sequences.entrySet()){
                if (entry1.getKey() <= entry2.getKey()){
                    continue;
                }
                int d = distance.apply(entry1.getValue(), entry2.getValue());
                if (d != -1 && d<=k){
                    result.add(new IntIntPair(entry1.getKey(), entry2.getKey()));
                }
            }
        }
        System.out.println("length = " +result.size());
        return result;
    }
}
