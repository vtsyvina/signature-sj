package by.bsu.model;

import java.util.Map;
import java.util.Queue;

/**
 * Structure to store sequence tree, where node key is prefix
 * Example:
 * sequences - AAAB, AAAC, AABA, CABA, CABA
 * one of tree options:
 *           ---------------root (key = '')---------------
 *           |                                           |
 *  -----node( key = 'AA')---------------            node(key = 'CABA', seq = 'CABA, CABA')
 *  |                                   |
 * node (key = 'AAA')----               node (key = 'AABA', seq = 'AABA')
 * |                     |
 * |                     node (key = 'AAAC', seq = 'AAAC')
 * |
 * node (key = 'AAAB', seq = 'AAAB')
 */
public class SequencesTree {
    public Node root;
    public static class Node{
        public String key;
        public Queue<Node> children;
        public Node parent;
        public Map<Integer, String> sequences;
    }
}
