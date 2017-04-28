package by.bsu.model;

import com.carrotsearch.hppc.LongArrayList;
import com.carrotsearch.hppc.ShortArrayList;

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
    public int l;
    public int maxChunks;
    public static class Node{
        public String key;
        public LongArrayList chunks;
        public Map<Long, ShortArrayList> grams;
        public Queue<Node> children;
        public Node parent;
        public Map<Integer, String> sequences;
    }
}
