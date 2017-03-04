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
        recursiveFillTree(0, sequences.size(), 0, tree.root, sequences);
        return tree;
    }

    private static void recursiveFillTree(int start, int end, int position,
                                          SequencesTree.Node currentNode, List<Pair<Integer, String>> sequences) {
        int left = start;
        while (left < end) {
            int max = 0;
            int maxBorder = -1;
            int maxI = 0;
            for (int i = 1; i < sequences.get(0).snd.length() - position; i++) {
                int border = start;
                while (border < end && sequences.get(border).snd.startsWith(sequences.get(left).snd.substring(0,position+i))){
                    border++;
                }
                if (i*(border-start+1) >= max){
                    max = i*(border-start+1);
                    maxBorder = border;
                    maxI = i;
                }
            }
            SequencesTree.Node node = new SequencesTree.Node();
            node.parent = currentNode;
            node.key = currentNode.key+sequences.get(left).snd.substring(position,position+maxI);
            if (maxBorder == -1){
                node.sequences = new HashMap<>();
                node.sequences.put(sequences.get(left).fst, sequences.get(left).snd);
            }else{
                node.children = new HashSet<>();
                recursiveFillTree(left, maxBorder, position+maxI, node, sequences);
            }
            left = maxBorder + 1;
        }
    }


    private static String greatestCommonPrefix(String a, String b) {
        int minLength = Math.min(a.length(), b.length());
        for (int i = 0; i < minLength; i++) {
            if (a.charAt(i) != b.charAt(i)) {
                return a.substring(0, i);
            }
        }
        return a.substring(0, minLength);
    }
}
