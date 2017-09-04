package by.bsu.util;

import by.bsu.model.IntStrPair;
import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;
import com.carrotsearch.hppc.LongArrayList;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Builds SequencesTree by given sample
 */
public class SequencesTreeBuilder {

    private static int maxChunks = 0;

    public static SequencesTree build(Sample sample) {
        SequencesTree tree = initEmptyTree();
        List<IntStrPair> sequences = new ArrayList<>();
        for (int i = 0; i < sample.sequences.length; i++) {
            sequences.add(new IntStrPair(i, sample.sequences[i]));
        }
        sequences.sort(Comparator.comparing(c -> c.r));
        recursiveFillTree(0, sequences.size() - 1, 0, tree.root, sequences);
        return tree;
    }

    public static SequencesTree build(Sample sample, int l) {
        SequencesTree tree = build(sample);
        maxChunks = 0;
        recursiveFillNodeChunks(tree.root, l);
        tree.l = l;
        tree.maxChunks = maxChunks;
        return tree;
    }

    private static SequencesTree initEmptyTree() {
        SequencesTree tree = new SequencesTree();
        tree.root = new SequencesTree.Node();
        tree.root.key = "";
        tree.root.parent = null;
        tree.root.children = new LinkedList<>();
        return tree;
    }

    public static void printTree(SequencesTree tree) {
        recursivePass(tree.root, 0);
    }

    private static void recursivePass(SequencesTree.Node node, int level) {
        if (node.children == null || node.children.isEmpty()) {
            StringBuilder str = new StringBuilder();
            for (int i = 0; i < level; i++) {
                str.append("   ");
            }
            System.out.println(str.toString() + node.key.length());
        } else {
            node.children.forEach(c -> recursivePass(c, level + 1));
            StringBuilder str = new StringBuilder();
            for (int i = 0; i < level; i++) {
                str.append("   ");
            }
            System.out.println(str.toString() + node.key.length() + " " + node.children.size());
        }
    }

    private static void recursiveFillNodeChunks(SequencesTree.Node node, int l) {
        node.chunks = new LongArrayList();
        node.grams = new HashMap<>();
        Utils.fillNodeGramsAndChunks(node, l);
        if (node.chunks.size() > maxChunks) {
            maxChunks = node.chunks.size();
        }
        if (node.children != null) {
            node.children.parallelStream().forEach(c -> recursiveFillNodeChunks(c, l));
        }
    }

    private static void recursiveFillTree(int start, int end, int prefixLength,
                                          SequencesTree.Node currentNode, List<IntStrPair> sequences) {
        int left = start;
        while (left <= end) {
            int maxLettersCount = 0;
            int maxBorder = -1;
            int maxPrefixLength = 0;
            for (int i = 1; i < sequences.get(0).r.length() - prefixLength + 1; i++) {
                int border = left;
                while (border < end && sequences.get(border + 1).r.startsWith(sequences.get(left).r.substring(0, prefixLength + i))) {
                    border++;
                }
                if (i * (border - left + 1) >= maxLettersCount) {
                    maxLettersCount = i * (border - left + 1);
                    maxBorder = border;
                    maxPrefixLength = i;
                }
            }
            SequencesTree.Node node = new SequencesTree.Node();
            node.parent = currentNode;
            node.key = currentNode.key + sequences.get(left).r.substring(prefixLength, prefixLength + maxPrefixLength);
            currentNode.children.add(node);
            if (prefixLength + maxPrefixLength == sequences.get(0).r.length()) {
                node.sequences = new HashMap<>();
                for (int i = left; i <= maxBorder; i++) {
                    node.sequences.put(sequences.get(i).l, sequences.get(i).r);
                }
            } else {
                node.children = new LinkedList<>();
                recursiveFillTree(left, maxBorder, prefixLength + maxPrefixLength, node, sequences);
            }
            left = maxBorder + 1;
        }
    }
}
