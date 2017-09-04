package by.bsu.util;

import by.bsu.model.KMerDict;
import by.bsu.model.Sample;
import com.carrotsearch.hppc.IntScatterSet;
import com.carrotsearch.hppc.LongScatterSet;

import java.util.HashMap;

import static by.bsu.util.Utils.convertLetterToDigit;

/**
 * Class to build KMerDict for given Sequence
 */
public class KMerDictBuilder {

    public static KMerDict getDict(Sample sample, int l){
        KMerDict result = new KMerDict();
        result.sampleName = sample.name;
        result.l = l;
        result.sequencesLength = sample.sequences[0].length();
        result.chunksCount = result.sequencesLength / l;

        result.sequencesNumber = sample.sequences.length;
        result.wholeSampleFixedPositionHashesList = new LongScatterSet[result.chunksCount];
        result.sequenceFixedPositionHashesList = new long[sample.sequences.length][];
        result.allHashesSet = new LongScatterSet();
        for (int i = 0; i < result.chunksCount; i++) {
            result.wholeSampleFixedPositionHashesList[i] = new LongScatterSet();
        }

        result.hashToSequencesMap = new HashMap<>(result.sequencesNumber);
        for (int seq = 0; seq < sample.sequences.length; seq++) {
            String sequence = sample.sequences[seq];
            long hashValue = 0;
            for (int j = 0; j < l; j++) {
                hashValue *= 4;
                hashValue += convertLetterToDigit(sequence.charAt(j));
            }
            result.sequenceFixedPositionHashesList[seq] = new long[result.chunksCount];
            result.sequenceFixedPositionHashesList[seq][0] = hashValue;
            result.wholeSampleFixedPositionHashesList[0].add(hashValue);
            if (!result.hashToSequencesMap.containsKey(hashValue)){
                result.hashToSequencesMap.put(hashValue, new IntScatterSet());
            }
            result.hashToSequencesMap.get(hashValue).add(seq);
            result.allHashesSet.add(hashValue);
            for (int j = 1; j < sequence.length() - l +1; j++) {
                hashValue -= ((long)convertLetterToDigit(sequence.charAt(j-1))) << 2 * (l-1);
                hashValue <<= 2;
                hashValue += convertLetterToDigit(sequence.charAt(j+l -1));
                if (!result.hashToSequencesMap.containsKey(hashValue)){
                    result.hashToSequencesMap.put(hashValue, new IntScatterSet());
                }
                result.hashToSequencesMap.get(hashValue).add(seq);
                result.allHashesSet.add(hashValue);
                if (j % l == 0){
                    result.sequenceFixedPositionHashesList[seq][j/l] = hashValue;
                    result.wholeSampleFixedPositionHashesList[j/l].add(hashValue);
                }
            }
        }
        return result;
    }

    public static KMerDict getDict(Sample sample, int l, int threshold){
        KMerDict result = getDict(sample, l);
        sample.consensus = Utils.consensus(sample.sequences);
        return result;
    }
}
