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
     *  how each position store map where key is l-chunk and value is array with marks of all sequences that contain this l-chunk
     *  each array has size of the sample. So chunksHashToSequencesMap[i].get(hash)[j] = 1 if j sequence has given hash as i chunk
     *  and 0 otherwise
     *  {
     *      0 -> {14323123 -> [0, 0, 1, 0, ..., 0, 1, 0,...], 36321322 -> [1, 0, 1, 1, 0, ... ,0, 1, ...]}
     *      1 -> {17481743 -> [0, 1, 1, 0, ..., 0, 1, 0,...], 96321422 -> [0, 0, 1, 1, 0, ... ,1, 0, ...]}
     *      ...
     *      fixedkMersCount -> {...}
     *  }
     *  
     */
    public Map<Long, int[]>[] chunksHashToSequencesMap;

    /**
     *  how each position store map where key is l-chunk and value is array of all sequences that contain this l-chunk
     *  each array in chunksHashToSequencesMapArray[i].get(hash) just stores list of all sequences indexes 
     *  that has given hash as i chunk
     *  {
     *      0 -> {14323123 -> [1, 4, 6, 7,...], 36321322 -> [1, 3, 6, 9,...]}
     *      1 -> {17481743 -> [2, 5, 6, 7,...], 96321422 -> [2, 5, 8, 9,...]}
     *      ...
     *      fixedkMersCount -> {...}
     *  }
     *
     */
    public Map<Long, int[]>[] chunksHashToSequencesMapArray;
}
