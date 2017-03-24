package by.bsu.algorithms;

import by.bsu.model.IntIntPair;
import by.bsu.model.Point;
import by.bsu.model.Points;
import by.bsu.model.Sample;

import by.bsu.util.HammingDistance;
import by.bsu.util.LevenshteinDistance;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.cursors.IntCursor;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Algorithm use points info to filter sequences from comparing
 * @deprecated Need to fix. Wrong output for some samples
 */
public class PointsMethod {

    public static Set<IntIntPair> run(Sample sample1, Sample sample2, Points points1, Points points2, int k){
        Set<IntIntPair> closePairs = new HashSet<>();
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        int comps = 0;
        boolean sameLength = sample1.sequences.values().iterator().next().length() ==
                sample2.sequences.values().iterator().next().length();
        for (Map.Entry<Point, IntSet> point1 : points1.pointSeqMap.entrySet()){
            for (Map.Entry<Point, IntSet> point2 : points2.pointSeqMap.entrySet()){
                if (lowerBoundEstimate(point1.getKey().arr, point2.getKey().arr) <= k){
                    for ( IntCursor s1 : point1.getValue()){
                        for (IntCursor s2 : point2.getValue()){
                            comps++;
                            if (sameLength && hammingDistance.apply(sample1.sequences.get(s1.value),
                                    sample2.sequences.get(s2.value)) <= k){
                                closePairs.add(new IntIntPair(s1.value, s2.value));
                                continue;
                            }
                            if (distance.apply(sample1.sequences.get(s1.value),
                                    sample2.sequences.get(s2.value)) != -1){
                                closePairs.add(new IntIntPair(s1.value, s2.value));
                            }
                        }
                    }
                }
            }
        }
        System.out.println("comps = "+comps);
        System.out.println("length = "+closePairs.size());
        return closePairs;
    }

    public static List<IntIntPair> run(Sample sample, Points points, int k){
        List<IntIntPair> closePairs = new ArrayList<>();
        LevenshteinDistance distance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        int comps = 0;
        Points toProcess = new Points(points);
        for (Map.Entry<Point, IntSet> point1 : points.pointSeqMap.entrySet()){
            for (Map.Entry<Point, IntSet> point2 : toProcess.pointSeqMap.entrySet()){
                if (lowerBoundEstimate(point1.getKey().arr, point2.getKey().arr) <= k){
                    for ( IntCursor s1 : point1.getValue()){
                        for (IntCursor s2 : point2.getValue()){
                            if (s1.value == s2.value){
                                continue;
                            }
                            comps++;
                            if (hammingDistance.apply(sample.sequences.get(s1.value),
                                    sample.sequences.get(s2.value)) <= k){
                                closePairs.add(new IntIntPair(s1.value, s2.value));
                                continue;
                            }
                            if (distance.apply(sample.sequences.get(s1.value),
                                    sample.sequences.get(s2.value)) != -1){
                                closePairs.add(new IntIntPair(s1.value, s2.value));
                            }
                        }
                    }
                }
            }
            toProcess.pointSeqMap.remove(point1.getKey());
        }
        System.out.println("comps = "+comps);
        System.out.println("length = "+closePairs.size());
        return closePairs;
    }

    public static int lowerBoundEstimate(int[] c1, int[] c2){
        int positive = 0;
        int negative = 0;
        for (int i = 0; i < 4; i++) {
            int dif = c1[i] - c2[i];
            if (dif > 0){
                positive += dif;
            } else{
                negative -= dif;
            }
        }
        return positive > negative ? positive : negative;
    }
}
