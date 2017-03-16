package by.bsu.algorithms;

import by.bsu.model.IntIntPair;
import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;
import by.bsu.util.HammingDistance;
import by.bsu.util.LevenshteinDistance;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by c5239200 on 3/5/17.
 */
public class TreeMethod {

    public static Set<IntIntPair> run(Sample sample, SequencesTree tree, int k) {
        Set<IntIntPair> result = ConcurrentHashMap.newKeySet();
        AtomicInteger count = new AtomicInteger(0);
        sample.sequences.entrySet().parallelStream().forEach(seq -> {
            recursiveDescent(seq, tree.root, k, result, count);
        });
        System.out.println("length = " + result.size());
        System.out.println("comps = " + count);
        return result;
    }

    public static Set<IntIntPair> runV2(Sample sample, SequencesTree tree, int k){
        Set<IntIntPair> result = ConcurrentHashMap.newKeySet();
        LevenshteinDistance levenshteinDistance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();

        tree.root.children.forEach(node -> recursiveDescentV2(node,
                tree.root.children,
                k, result, levenshteinDistance, hammingDistance));
        return result;
    }

    private static void recursiveDescent(Map.Entry<Integer, String> entry, SequencesTree.Node node, int k, Set<IntIntPair> result, AtomicInteger count) {
        LevenshteinDistance levenshteinDistance = new LevenshteinDistance(k);
        HammingDistance hammingDistance = new HammingDistance();
        count.incrementAndGet();
        if (hammingDistance.apply(entry.getValue().substring(0, node.key.length()), node.key) <= k
                || levenshteinDistance.apply(entry.getValue().substring(0, node.key.length()), node.key) != -1) {
            if (node.children != null && !node.children.isEmpty()) {
                node.children.forEach(n -> recursiveDescent(entry, n, k, result, count));
            }
            if (node.sequences != null) {
                node.sequences.entrySet().forEach(s ->
                {
                    if (!(entry.getKey() <= s.getKey())) {
                        if (hammingDistance.apply(entry.getValue(), s.getValue()) <= k) {
                            result.add(new IntIntPair(entry.getKey(), s.getKey()));
                        } else {
                            //count.incrementAndGet();
                            if (levenshteinDistance.apply(entry.getValue(), s.getValue()) != -1) {
                                result.add(new IntIntPair(entry.getKey(), s.getKey()));
                            }
                        }
                    }
                });
            }
        }
    }

    private static void recursiveDescentV2(SequencesTree.Node node, Set<SequencesTree.Node> toCheck, int k,
                                           Set<IntIntPair> result, LevenshteinDistance levenshtein, HammingDistance hamming){
        if (node.sequences == null){
            Set<SequencesTree.Node> newToCheck = new HashSet<>();
            toCheck.forEach( check -> newToCheck.addAll(calculateToCheckForGivenNode(node, check, k , levenshtein, hamming)));
            if (node.children != null){
                node.children.forEach( child -> recursiveDescentV2(child, newToCheck, k, result, levenshtein, hamming));
            }
        } else {
            Queue<SequencesTree.Node> nodesToVisit = new LinkedList<>(toCheck);
            while (!nodesToVisit.isEmpty()){
                SequencesTree.Node currentNode = nodesToVisit.poll();
                if (currentNode.children != null){
                    int min = currentNode.key.length() < node.key.length() ? currentNode.key.length() : node.key.length();
                    int h = hamming.apply(currentNode.key.substring(0, min), node.key.substring(0, min));
                    int l = levenshtein.apply(currentNode.key.substring(0, min), node.key.substring(0, min));
                    if ( h <= k ||  l != -1) {
                        nodesToVisit.addAll(currentNode.children);
                    }
                } else if (currentNode != node){
                    int h = hamming.apply(currentNode.key, node.key);
                    int l = levenshtein.apply(currentNode.key, node.key);
                    if ( h <= k ||  l != -1){
                        node.sequences.keySet().forEach( index ->
                                currentNode.sequences.keySet().forEach(
                                        index2 -> result.add(new IntIntPair(index, index2))));
                    }
                } else if(currentNode.sequences.size() > 1){
                    node.sequences.keySet().forEach( index ->
                            currentNode.sequences.keySet().forEach(
                                    index2 -> {
                                        if (index != index2)
                                            result.add(new IntIntPair(index, index2));
                                    }));

                }
            }
        }
        //node.parent.children.remove(node);
    }

    /**
     * Calculates set of nodes to check for currentNode for toCheck node and its children
     * return set ot nodes that satisfies  condition:
     * currentNode.key.length() < toCheck.key.length() && levenshtein(currentNode.key, toCheck.key[0:currentNode.key.length])<= k)
     * @param currentNode given node
     * @param toCheck current node to check
     * @param k threshold
     * @return set of nodes with key length bigger than currentNode.key and with distance less than k
     */
    private static Set<SequencesTree.Node> calculateToCheckForGivenNode(SequencesTree.Node currentNode, SequencesTree.Node toCheck, int k, LevenshteinDistance levenshtein, HammingDistance hamming){
        Set<SequencesTree.Node> result = new HashSet<>();
        int min = currentNode.key.length() < toCheck.key.length() ? currentNode.key.length() : toCheck.key.length();
        int h = hamming.apply(currentNode.key.substring(0, min), toCheck.key.substring(0, min));
        int l = levenshtein.apply(currentNode.key.substring(0, min), toCheck.key.substring(0, min));
        if ( h <= k ||  l != -1){
            if (currentNode.key.length() <= toCheck.key.length()){
                result.add(toCheck);
            } else{
                if (toCheck.children != null){
                    toCheck.children.forEach( child -> result.addAll(calculateToCheckForGivenNode(currentNode, child, k ,levenshtein, hamming)));
                }
            }
        }
        return result;
    }
}
