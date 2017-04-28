package by.bsu.algorithms;

import by.bsu.model.IntIntPair;
import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;
import by.bsu.start.Start;
import by.bsu.util.HammingDistance;
import by.bsu.util.LevenshteinDistance;
import com.carrotsearch.hppc.ShortArrayList;
import com.carrotsearch.hppc.cursors.ShortCursor;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import static by.bsu.util.Utils.numbers;

/**
 * Created by c5239200 on 3/5/17.
 */
public class TreeMethod {

    //debug variable for all comparisons
    private static int allComps = 0;
    //debug variable to count levenshtein comparisons
    private static int levenshteinSeqComp = 0;
    //debug variable to now how many cpmparisons were reduces by singnature check
    private static int signatureReduce = 0;
    //store result length
    private static int resultLength = 0;
    //service variable for periodic write to file
    private static int previousWrite = 0;
    // builder that stores intermediate result before writing to file
    private static StringBuilder builder;
    // algorithm output
    private static Path path;
    //global algorithm variables
    private static int l = 0;
    private static int maxChunks = 0;
    private static int k = 0;
    private static LevenshteinDistance levenshtein;
    private static HammingDistance hamming;

    public static Set<IntIntPair> run(Sample sample, SequencesTree tree, int k) {

        Set<IntIntPair> result = ConcurrentHashMap.newKeySet();
        AtomicInteger count = new AtomicInteger(0);
        sample.sequences.entrySet().parallelStream().forEach(seq -> {
            recursiveDescent(seq, tree.root, k, result, count);
        });
        System.out.println("length = " + result.size());
        System.out.println("allComps = " + count);
        return result;
    }

    /**
     * Only for samples with the same length
     */
    public static long runV2(Sample sample, SequencesTree tree, int k) throws IOException {
        System.out.println("Start Tree method for " + sample.name + " k=" + k);
        allComps = 0;
        levenshteinSeqComp = 0;
        signatureReduce = 0;
        l = tree.l;
        resultLength = 0;
        previousWrite = 0;
        TreeMethod.k = k;
        maxChunks = tree.maxChunks;
        builder = new StringBuilder();
        path = Start.getOutputFilename(sample, "tree");

        levenshtein = new LevenshteinDistance(k);
        hamming = new HammingDistance();
        while (!tree.root.children.isEmpty()) {
            recursiveDescentV2(tree.root.children.peek(),
                    tree.root.children,
                    k);
        }
        System.out.println();
        System.out.println("comps = " + (allComps + levenshteinSeqComp));
        System.out.println("travel allComps = " + (allComps));
        System.out.println("levenshtein = " + (levenshteinSeqComp));
        System.out.println("signatureReduce = " + signatureReduce);
        System.out.println("length = " + resultLength);
        return resultLength;
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

    private static void recursiveDescentV2(SequencesTree.Node node, Queue<SequencesTree.Node> toCheck, int k) throws IOException {
        if (node.sequences == null) {
            Queue<SequencesTree.Node> newToCheck = new LinkedList<>();
            toCheck.forEach(check -> newToCheck.addAll(calculateToCheckForGivenNode(node, check, levenshtein, hamming)));
            if (node.children != null) {
                while (!node.children.isEmpty()) {
                    recursiveDescentV2(node.children.peek(), newToCheck, k);
                }

            }
        } else {
            Queue<SequencesTree.Node> nodesToVisit = new LinkedList<>(toCheck);
            while (!nodesToVisit.isEmpty()) {
                SequencesTree.Node currentNode = nodesToVisit.poll();
                if (currentNode.children != null) {
                    walkLevelDeeper(node, k, nodesToVisit, currentNode);
                } else if (currentNode != node) {
                    compareSequencesFromDifferentNodes(node, k, currentNode);
                } else if (currentNode.sequences.size() > 1) {
                    addSequencesFromSameNode(node, currentNode);

                }
            }
        }
        node.parent.children.remove(node);
    }

    /**
     * If node is in toCheckList but doesn't contain sequences we need to go further to get sequences from its children
     */
    private static void walkLevelDeeper(SequencesTree.Node node, int k, Queue<SequencesTree.Node> nodesToVisit, SequencesTree.Node currentNode) {
        int min = currentNode.key.length() < node.key.length() ? currentNode.key.length() : node.key.length();
        allComps++;
        if (hamming.apply(currentNode.key.substring(0, min), node.key.substring(0, min)) <= k
                || signatureCheck(currentNode, node)
                || levenshtein.apply(currentNode.key.substring(0, min), node.key.substring(0, min)) != -1) {
            nodesToVisit.addAll(currentNode.children);
        }
    }

    /**
     * Simply compares all sequences from first node with all sequences from second node
     */
    private static void compareSequencesFromDifferentNodes(SequencesTree.Node node, int k, SequencesTree.Node currentNode) throws IOException {
        if (hamming.apply(currentNode.key, node.key) <= k) {
            node.sequences.keySet().forEach(index ->
                    currentNode.sequences.keySet().forEach(
                            index2 -> {
                                resultLength++;
                                builder.append(numbers.get(index)).append(" ").append(numbers.get(index2)).append("\n");
                            }
                    )
            );
        } else {
            if (!signatureCheck(currentNode, node)) {
                return;
            }
            levenshteinSeqComp++;
            if (levenshtein.apply(currentNode.key, node.key) != -1) {
                node.sequences.keySet().forEach(index ->
                        currentNode.sequences.keySet().forEach(
                                index2 -> {
                                    resultLength++;
                                    builder.append(numbers.get(index)).append(" ").append(numbers.get(index2)).append("\n");
                                }
                        )
                );
            }
        }
        appendResultToFile();
    }

    /**
     * If sample has several equal sequences we need to add all their pair combinations
     */
    private static void addSequencesFromSameNode(SequencesTree.Node node, SequencesTree.Node currentNode) throws IOException {
        node.sequences.keySet().forEach(index ->
                currentNode.sequences.keySet().forEach(
                        index2 -> {
                            if (index < index2) {
                                resultLength++;
                                builder.append(numbers.get(index)).append(" ").append(numbers.get(index2)).append("\n");
                            }

                        }
                )
        );
        appendResultToFile();
    }

    /**
     * Calculates set of nodes to check for currentNode for toCheck node and its children
     * return set ot nodes that satisfies  condition:
     * currentNode.key.length() < toCheck.key.length() && levenshtein(currentNode.key, toCheck.key[0:currentNode.key.length])<= k)
     *
     * @param currentNode given node
     * @param toCheck     current node to check
     * @return set of nodes with key length bigger than currentNode.key and with distance less than k
     */
    private static Set<SequencesTree.Node> calculateToCheckForGivenNode(SequencesTree.Node currentNode, SequencesTree.Node toCheck, LevenshteinDistance levenshtein, HammingDistance hamming) {
        Set<SequencesTree.Node> result = new HashSet<>();
        int min = currentNode.key.length() < toCheck.key.length() ? currentNode.key.length() : toCheck.key.length();
        allComps++;
        if (hamming.apply(currentNode.key.substring(0, min), toCheck.key.substring(0, min)) <= k
                || signatureCheck(currentNode, toCheck)
                || levenshtein.apply(currentNode.key.substring(0, min), toCheck.key.substring(0, min)) != -1) {
            if (currentNode.key.length() <= toCheck.key.length()) {
                result.add(toCheck);
            } else {
                if (toCheck.children != null) {
                    toCheck.children.forEach(child -> result.addAll(calculateToCheckForGivenNode(currentNode, child, levenshtein, hamming)));
                }
            }
        }
        return result;
    }

    /**
     * Checks if one of node contains enough grams for chunks of second node
     */
    private static boolean signatureCheck(SequencesTree.Node n1, SequencesTree.Node n2) {
        SequencesTree.Node min = n1.key.length() < n2.key.length() ? n1 : n2;
        SequencesTree.Node max = n1.key.length() < n2.key.length() ? n2 : n1;
        int missMatches = 0;
        for (int i = 0; i < min.chunks.size(); i++) {
            if (!fitTheThreshold(max.grams.get(min.chunks.get(i)), i * l, k)) {
                missMatches++;
            }
        }
        signatureReduce++;
        return missMatches <= maxChunks - k;
    }

    private static boolean fitTheThreshold(ShortArrayList list, int position, int threshold) {
        if (list == null) {
            return false;
        }
        for (ShortCursor shortCursor : list) {
            if (shortCursor.value >= position - threshold && shortCursor.value <= position + threshold) {
                return true;
            }
        }
        return false;
    }

    private static void appendResultToFile() throws IOException {
        if (resultLength - previousWrite > 5000) {
            Files.write(path, builder.toString().getBytes(), StandardOpenOption.APPEND);
            builder = new StringBuilder();
            previousWrite = resultLength;
            System.out.print("\r" + resultLength);
        }
    }
}