package by.bsu.model;

import com.carrotsearch.hppc.*;

import java.util.HashMap;
import java.util.Map;

/**
 * Data container to store hashes of samples k-mers
 */
public class KMerDict {
    /**
     *  how each possible l-mer in sample store Set of sequences that contain given l-mer
     *  14323123 -> set(1, 4, 6, 7,...)
     *  36321322 -> set(1, 3, 6, 9,...)
     */
    public Map<Long, IntSet> hashToSequencesMap;

    /**
     * Contains only all possible l-mers in sample
     */
    public LongSet allHashesSet;

    /**
     * for each sequence from sample stores hashes for l-mers on 0, 1*l, 2*l,... positions(in the amount of fixedkMersCount)
     */
    public Map<Integer, long[]> sequenceFixedPositionHashesList;
    /**
     * for each position 0, 1*l, 2*l,... store set of hashes that occur at any sequence from sample
     */
    public LongSet[] wholeSampleFixedPositionHashesList;

    /**
     * Store distance from each sequence to consensus string
     */
    public Map<Integer, Integer> consensusDistances;

    /**
     * Length of l-mers partition
     */
    public int l;
    /**
     * Count of l-mers in partition
     */
    public int fixedkMersCount;
    /**
     * Length of taken separately sequence
     */
    public int sequencesLength;
    /**
     * The number of sequences in sample
     */
    public int sequencesNumber;
    /**
     * Sample's name
     */
    public String sampleName;
}
