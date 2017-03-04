package by.bsu.model;

import java.util.Map;
import java.util.Set;

/**
 * Created by c5239200 on 2/28/17.
 */
public class SequencesTree {
    public Node root;
    public static class Node{
        public String key;
        public Set<Node> children;
        public Node parent;
        public Map<Integer, String> sequences;
    }
}
