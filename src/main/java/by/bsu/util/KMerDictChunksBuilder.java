package by.bsu.util;

import static by.bsu.util.Utils.getHashValue;

import java.util.HashMap;
import java.util.Map;

import com.carrotsearch.hppc.IntScatterSet;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongScatterSet;
import com.carrotsearch.hppc.cursors.IntCursor;

import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;

/**
 * Created by c5239200 on 6/26/17.
 */
public class KMerDictChunksBuilder {

    public static KMerDictChunks getDict(Sample sample, int l){
        KMerDictChunks result = new KMerDictChunks();
        result.sampleName = sample.name;
        result.l = l;
        result.sequencesLength = sample.sequences[0].length();
        result.fixedkMersCount = result.sequencesLength / l;

        result.sequencesNumber = sample.sequences.length;
        result.wholeSampleFixedPositionHashesList = new LongScatterSet[result.fixedkMersCount];
        result.sequenceFixedPositionHashesList = new long[sample.sequences.length][];
        result.allHashesSet = new LongScatterSet();
        for (int i = 0; i < result.fixedkMersCount; i++) {
            result.wholeSampleFixedPositionHashesList[i] = new LongScatterSet();
        }

        result.chunksHashToSequencesMap = new HashMap[result.fixedkMersCount];
        result.chunksHashToSequencesMapArray = new HashMap[result.fixedkMersCount];
        for (int i = 0; i < result.fixedkMersCount; i++) {
            result.chunksHashToSequencesMap[i] = new HashMap<>();
            result.chunksHashToSequencesMapArray[i] = new HashMap<>();
        }
        for (int seq = 0; seq < sample.sequences.length; seq++) {
            result.sequenceFixedPositionHashesList[seq] = new long[result.fixedkMersCount];
            for (int i = 0; i < result.fixedkMersCount; i++) {
                long hashValue = getHashValue(i*l, l, sample.sequences[seq]);
                result.sequenceFixedPositionHashesList[seq][i] = hashValue;
                result.wholeSampleFixedPositionHashesList[i].add(hashValue);
                if (!result.chunksHashToSequencesMap[i].containsKey(hashValue)){
                    result.chunksHashToSequencesMap[i].put(hashValue, new int[sample.sequences.length]);
                }
                result.chunksHashToSequencesMap[i].get(hashValue)[seq] = 1;
                result.allHashesSet.add(hashValue);
            }
        }
        for (int i = 0; i < result.chunksHashToSequencesMap.length; i++) {
            for (Map.Entry<Long, int[]> entry : result.chunksHashToSequencesMap[i].entrySet()) {
                int[] tmp = new int[entry.getValue().length];
                int j = 0;
                for (int k = 0; k < entry.getValue().length; k++) {
                    if (entry.getValue()[k] == 1){
                        tmp[j++] = k; 
                    }
                }
                result.chunksHashToSequencesMapArray[i].put(entry.getKey(), tmp);
            }
        }
        return result;
    }
}
