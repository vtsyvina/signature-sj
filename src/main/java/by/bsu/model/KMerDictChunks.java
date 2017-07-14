package by.bsu.model;

import java.util.Map;

import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongSet;

/**
 * Created by c5239200 on 6/26/17.
 * Another model to store sample data that operates only with chunks
 */
public class KMerDictChunks extends AbstractKMerDict{
    
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
}
