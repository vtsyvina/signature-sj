package by.bsu.util;

import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;
import com.sun.tools.javac.util.Pair;

import java.util.*;

/**
 * Created by c5239200 on 2/28/17.
 */
public class SequencesTreeBuilder {

    public static SequencesTree build(Sample sample) {
        SequencesTree tree = new SequencesTree();
        tree.root = new SequencesTree.Node();
        tree.root.key = "";
        tree.root.parent = null;
        tree.root.children = new HashSet<>();

        List<Pair<Integer, String>> sequences = new ArrayList<>();
        sample.sequences.entrySet().forEach(s -> sequences.add(new Pair<>(s.getKey(), s.getValue())));
        sequences.sort(Comparator.comparing(c -> c.snd));
        recursiveFillTree(0, sequences.size()-1, 0, tree.root, sequences);
        return tree;
    }

    private static void recursiveFillTree(int start, int end, int position,
                                          SequencesTree.Node currentNode, List<Pair<Integer, String>> sequences) {
        int left = start;
        while (left <= end) {
            int max = 0;
            int maxBorder = -1;
            int maxI = 0;
            for (int i = 1; i < sequences.get(0).snd.length() - position +1; i++) {
                int border = left;
                while (border < end && sequences.get(border+1).snd.startsWith(sequences.get(left).snd.substring(0,position+i))){
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
            node.key = currentNode.key+sequences.get(left).snd.substring(position,position+maxI);
            currentNode.children.add(node);
            if (position+maxI == sequences.get(0).snd.length()){
                node.sequences = new HashMap<>();
                for (int i = left; i <= maxBorder; i++) {
                    node.sequences.put(sequences.get(i).fst, sequences.get(i).snd);
                }
            }else{
                node.children = new HashSet<>();
                recursiveFillTree(left, maxBorder, position+maxI, node, sequences);
            }
            left = maxBorder + 1;
        }
    }
}
