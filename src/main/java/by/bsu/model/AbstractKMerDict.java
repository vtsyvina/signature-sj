package by.bsu.model;

import java.util.Map;

import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongSet;

/**
 * Created by c5239200 on 7/13/17.
 */
public abstract class AbstractKMerDict {

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
     * for each sequence from sample stores hashes for l-mers on 0, 1*l, 2*l,... positions(in the amount of chunksCount)
     */
    public long[][] sequenceFixedPositionHashesList;
    /**
     * for each position 0, 1*l, 2*l,... store set of hashes that occur at any sequence from sample
     */
    public LongSet[] wholeSampleFixedPositionHashesList;

    /**
     * Length of l-mers partition
     */
    public int l;
    /**
     * Count of l-mers in partition
     */
    public int chunksCount;
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
