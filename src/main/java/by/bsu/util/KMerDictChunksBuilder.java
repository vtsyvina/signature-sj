package by.bsu.util;

import static by.bsu.util.Utils.convertLetterToDigit;

import java.util.HashMap;
import java.util.Map;

import com.carrotsearch.hppc.IntScatterSet;
import com.carrotsearch.hppc.LongScatterSet;

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
        result.sequencesLength = sample.sequences.values().iterator().next().length();
        result.fixedkMersCount = result.sequencesLength / l;

        result.sequencesNumber = sample.sequences.size();
        result.wholeSampleFixedPositionHashesList = new LongScatterSet[result.fixedkMersCount];
        result.sequenceFixedPositionHashesList = new HashMap<>();
        result.allHashesSet = new LongScatterSet();
        for (int i = 0; i < result.fixedkMersCount; i++) {
            result.wholeSampleFixedPositionHashesList[i] = new LongScatterSet();
        }

        result.chunksHashToSequencesMap = new HashMap[result.fixedkMersCount];
        for (int i = 0; i < result.fixedkMersCount; i++) {
            result.chunksHashToSequencesMap[i] = new HashMap<>();
        }
        for (Map.Entry<Integer, String> entry : sample.sequences.entrySet()) {
            result.sequenceFixedPositionHashesList.put(entry.getKey(), new long[result.fixedkMersCount]);
            for (int i = 0; i < result.fixedkMersCount; i++) {
                long hashValue = getHashValue(i*l, l, entry);
                result.sequenceFixedPositionHashesList.get(entry.getKey())[i] = hashValue;
                result.wholeSampleFixedPositionHashesList[i].add(hashValue);
                if (!result.chunksHashToSequencesMap[i].containsKey(hashValue)){
                    result.chunksHashToSequencesMap[i].put(hashValue, new IntScatterSet());
                }
                result.chunksHashToSequencesMap[i].get(hashValue).add(entry.getKey());
                result.allHashesSet.add(hashValue);
            }
        }
        return result;
    }

    private static long getHashValue(int position, int l, Map.Entry<Integer, String> entry) {
        long hashValue = 0;
        for (int j = 0; j < l; j++) {
            hashValue *= 4;
            hashValue += convertLetterToDigit(entry.getValue().charAt(position+j));
        }
        return hashValue;
    }
}
