package by.bsu.model;

import java.util.Map;

import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongSet;

/**
 * Created by c5239200 on 6/26/17.
 * Another model to store sample data that operates only with chunks
 */
public class KMerDictChunks {
    
    /**
     *  how each position store map where key is l-chunk and value is set of all sequences that contain this l-chunk
     *  {
     *      0 -> {14323123 -> set(1, 4, 6, 7,...), 36321322 -> set(1, 3, 6, 9,...)}
     *      1 -> {17481743 -> set(2, 5, 6, 7,...), 96321422 -> set(2, 5, 8, 9,...)}
     *      ...
     *      fixedkMersCount -> {...}
     *  }
     *  
     */
    public Map<Long, IntSet>[] chunksHashToSequencesMap;

    /**
     * Contains only all possible l-chunks in sample
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
