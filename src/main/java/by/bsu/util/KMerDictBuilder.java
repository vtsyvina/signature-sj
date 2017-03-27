package by.bsu.util;

import by.bsu.model.KMerDict;
import by.bsu.model.Sample;
import com.carrotsearch.hppc.*;

import java.util.HashMap;
import java.util.Map;

import static by.bsu.util.Utils.convertLetterToDigit;

/**
 * Class to build KMerDict for given Sequence
 */
public class KMerDictBuilder {

    public static KMerDict getDict(Sample sample, int l){
        KMerDict result = new KMerDict();
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

        result.hashToSequencesMap = new HashMap<>(result.sequencesNumber);
        for (Map.Entry<Integer, String> entry : sample.sequences.entrySet()) {
            long hashValue = 0;
            for (int j = 0; j < l; j++) {
                hashValue *= 4;
                hashValue += convertLetterToDigit(entry.getValue().charAt(j));
            }
            result.sequenceFixedPositionHashesList.put(entry.getKey(), new long[result.fixedkMersCount]);
            result.sequenceFixedPositionHashesList.get(entry.getKey())[0] = hashValue;
            result.wholeSampleFixedPositionHashesList[0].add(hashValue);
            if (!result.hashToSequencesMap.containsKey(hashValue)){
                result.hashToSequencesMap.put(hashValue, new IntScatterSet());
            }
            result.hashToSequencesMap.get(hashValue).add(entry.getKey());
            result.allHashesSet.add(hashValue);
            for (int j = 1; j < entry.getValue().length() - l +1; j++) {
                hashValue -= ((long)convertLetterToDigit(entry.getValue().charAt(j-1))) << 2 * (l-1);
                hashValue <<= 2;
                hashValue += convertLetterToDigit(entry.getValue().charAt(j+l -1));
                if (!result.hashToSequencesMap.containsKey(hashValue)){
                    result.hashToSequencesMap.put(hashValue, new IntScatterSet());
                }
                result.hashToSequencesMap.get(hashValue).add(entry.getKey());
                result.allHashesSet.add(hashValue);
                if (j % l == 0){
                    result.sequenceFixedPositionHashesList.get(entry.getKey())[j/l] = hashValue;
                    result.wholeSampleFixedPositionHashesList[j/l].add(hashValue);
                }
            }
        }
        return result;
    }


}
