package by.bsu.util;

import by.bsu.model.IntStrPair;
import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;

import java.util.*;

/**
 * Builds SequencesTree by given sample
 */
public class SequencesTreeBuilder {

    public static SequencesTree build(Sample sample) {
        SequencesTree tree = new SequencesTree();
        tree.root = new SequencesTree.Node();
        tree.root.key = "";
        tree.root.parent = null;
        tree.root.children = new LinkedList<>();

        List<IntStrPair> sequences = new ArrayList<>();
        sample.sequences.entrySet().forEach(s -> sequences.add(new IntStrPair(s.getKey(), s.getValue())));
        sequences.sort(Comparator.comparing(c -> c.r));
        recursiveFillTree(0, sequences.size()-1, 0, tree.root, sequences);
        return tree;
    }

    private static void recursiveFillTree(int start, int end, int position,
                                          SequencesTree.Node currentNode, List<IntStrPair> sequences) {
        int left = start;
        while (left <= end) {
            int max = 0;
            int maxBorder = -1;
            int maxI = 0;
            for (int i = 1; i < sequences.get(0).r.length() - position +1; i++) {
                int border = left;
                while (border < end && sequences.get(border+1).r.startsWith(sequences.get(left).r.substring(0,position+i))){
                    border++;
                }
                if (i*(border-left+1) >= max){
                    max = i*(border-left+1);
                    maxBorder = border;
                    maxI = i;
                }
            }
            SequencesTree.Node node = new SequencesTree.Node();
            node.parent = currentNode;
            node.key = currentNode.key+sequences.get(left).r.substring(position,position+maxI);
            currentNode.children.add(node);
            if (position+maxI == sequences.get(0).r.length()){
                node.sequences = new HashMap<>();
                for (int i = left; i <= maxBorder; i++) {
                    node.sequences.put(sequences.get(i).l, sequences.get(i).r);
                }
            }else{
                node.children = new LinkedList<>();
                recursiveFillTree(left, maxBorder, position+maxI, node, sequences);
            }
            left = maxBorder + 1;
        }
    }
}
