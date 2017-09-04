package by.bsu.util;

import static by.bsu.util.Utils.getHashValue;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;

import com.carrotsearch.hppc.LongScatterSet;

import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;

/**
 * Created by c5239200 on 6/26/17.
 */
public class KMerDictChunksBuilder {

    private static final double eps = 0.000_001;

    /**
     * Calculates equal size chunks
     */
    public static KMerDictChunks getDict(Sample sample, int l) {
        KMerDictChunks result = new KMerDictChunks();
        result.l = l;
        result.sequencesLength = sample.sequences[0].length();
        result.chunksCount = result.sequencesLength / l;
        return get(result, sample, r -> r * l, r -> l);
    }

    /**
     * Calculates chunks based on sample profile entropy
     */
    public static KMerDictChunks getDict(Sample sample, int chunksCount, double[][] profile) {
        KMerDictChunks result = new KMerDictChunks();
        result.chunksCount = chunksCount;
        int[] chunkEnds = getChunksEnds(chunksCount, sample.sequences[0].length(), profile);
        return get(result, sample, r -> r == 0 ? 0 : chunkEnds[r - 1], r -> r == 0 ? chunkEnds[r] : chunkEnds[r] - chunkEnds[r - 1]);
    }

    private static KMerDictChunks get(KMerDictChunks result, Sample sample, Function<Integer, Integer> leftSymbol, Function<Integer, Integer> rightSymbol) {
        result.sampleName = sample.name;
        result.sequencesNumber = sample.sequences.length;
        result.wholeSampleFixedPositionHashesList = new LongScatterSet[result.chunksCount];
        result.sequenceFixedPositionHashesList = new long[sample.sequences.length][];
        result.allHashesSet = new LongScatterSet();
        for (int i = 0; i < result.chunksCount; i++) {
            result.wholeSampleFixedPositionHashesList[i] = new LongScatterSet();
        }

        result.chunksHashToSequencesMap = new HashMap[result.chunksCount];
        result.chunksHashToSequencesMapArray = new HashMap[result.chunksCount];
        for (int i = 0; i < result.chunksCount; i++) {
            result.chunksHashToSequencesMap[i] = new HashMap<>();
            result.chunksHashToSequencesMapArray[i] = new HashMap<>();
        }
        for (int seq = 0; seq < sample.sequences.length; seq++) {
            result.sequenceFixedPositionHashesList[seq] = new long[result.chunksCount];
            for (int i = 0; i < result.chunksCount; i++) {
                long hashValue = getHashValue(leftSymbol.apply(i), rightSymbol.apply(i), sample.sequences[seq]);
                result.sequenceFixedPositionHashesList[seq][i] = hashValue;
                result.wholeSampleFixedPositionHashesList[i].add(hashValue);
                if (!result.chunksHashToSequencesMap[i].containsKey(hashValue)) {
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
                    if (entry.getValue()[k] == 1) {
                        tmp[j++] = k;
                    }
                }
                result.chunksHashToSequencesMapArray[i].put(entry.getKey(), tmp);
            }
        }
        return result;
    }

    private static int[] getChunksEnds(int chunksCount, int sequencesLength, double[][] profile) {
        double[] entropy = new double[sequencesLength];
        //calculate profile entropy for each position
        for (int i = 0; i < sequencesLength; i++) {
            for (int j = 0; j < 4; j++) {
                if (profile[j][i] < eps) {
                    continue;
                }
                entropy[i] -= profile[j][i] * Math.log(profile[j][i]);
            }
        }
        double entSum = Arrays.stream(entropy).sum();
        int[] chunksEnds = new int[chunksCount];
        int currentPosition = 0;
        int prevPosition = 0;
        double currentSum = 0;
        double entropyPerChunk = entSum / chunksCount;
        // calculate chunks such that each have equal sum of entropy unless it has length of 32
        for (int i = 0; i < chunksCount; i++) {
            if (i == chunksCount - 1) {
                chunksEnds[i] = sequencesLength - 1;
                continue;
            }
            while ((Math.abs(currentSum - entropyPerChunk) > Math.abs(currentSum - entropyPerChunk + entropy[currentPosition]) || entropy[currentPosition] < eps)
                    && currentPosition - prevPosition < 32) {
                currentSum += entropy[currentPosition];
                currentPosition++;
            }
            chunksEnds[i] = currentPosition - 1;
            //recalculate entropy per chunk because current chunk has fewer entropy sum than it should
            if ((i == 0 && chunksEnds[i] == 31) || (i != 0 && chunksEnds[i] - chunksEnds[i - 1] == 32)) {
                entropyPerChunk += (entropyPerChunk - currentSum) / (chunksCount - i - 1);
            }
            currentSum = 0;
            prevPosition = currentPosition;
        }
        return chunksEnds;
    }
}
