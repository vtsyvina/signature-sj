package by.bsu.algorithms;

import by.bsu.model.Pair;
import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;
import by.bsu.util.HammingDistance;
import by.bsu.util.LevenshteinDistance;

import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by c5239200 on 3/5/17.
 */
public class TreeMethod {

    public static Set<Pair> run(Sample sample, SequencesTree tree, int k) {
        Set<Pair> result = ConcurrentHashMap.newKeySet();
        AtomicInteger count = new AtomicInteger(0);
        AtomicInteger reduce = new AtomicInteger(0);
        sample.sequences.entrySet().parallelStream().forEach(seq -> {
            recursiveDescent(seq, tree.root, k, result, count, reduce);
        });
        System.out.println("length = " + result.size());
        System.out.println("comps = " + count);
        System.out.println("reduce = " + reduce);
        return result;
    }

    private static void recursiveDescent(Map.Entry<Integer, String> entry, SequencesTree.Node node, int k, Set<Pair> result, AtomicInteger count, AtomicInteger reduce) {
        LevenshteinDistance levenshteinDistance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        count.incrementAndGet();
        if (hammingDistance.apply(entry.getValue().substring(0, node.key.length()), node.key) <= k
                || levenshteinDistance.apply(entry.getValue().substring(0, node.key.length()), node.key) != -1) {
            if (node.children != null && !node.children.isEmpty()) {
                node.children.forEach(n -> recursiveDescent(entry, n, k, result, count, reduce));
            }
            if (node.sequences != null) {
                node.sequences.entrySet().forEach(s ->
                {
                    if (!s.getKey().equals(entry.getKey())) {
                        if (hammingDistance.apply(entry.getValue(), s.getValue()) <= k) {
                            result.add(new Pair(entry.getKey(), s.getKey()));
                        } else {
                            //count.incrementAndGet();
                            if (levenshteinDistance.apply(entry.getValue(), s.getValue()) != -1) {
                                result.add(new Pair(entry.getKey(), s.getKey()));
                            }
                        }
                    }
                });
            }
        } else {
            //c(node, reduce);
        }
    }

    private static void c(SequencesTree.Node node, AtomicInteger reduce){
        if (node.children != null && !node.children.isEmpty()) {
            node.children.forEach(n -> c(n, reduce));
        }
        if (node.sequences != null) {
            reduce.addAndGet(node.sequences.size());
        }
    }
}
