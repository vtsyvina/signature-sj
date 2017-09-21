package by.bsu.util;

import com.carrotsearch.hppc.LongSet;
import com.carrotsearch.hppc.cursors.LongCursor;

import by.bsu.model.AbstractKMerDict;

/**
 * Created by c5239200 on 7/13/17.
 */
public class AlgorithmUtils {

    /**
     * Calculates how many l-mers from fixed positions from first dictionary are in second dictionary
     */
    public static int calculateCoincidences(AbstractKMerDict dict1, AbstractKMerDict dict2) {
        int kMerCoincidences = 0;
        for (LongSet positionHashes : dict1.wholeSampleChunksHashesList) {
            for (LongCursor hash : positionHashes) {
                if (dict2.allHashesSet.contains(hash.value)) {
                    kMerCoincidences++;
                    break;
                }
            }
        }
        return kMerCoincidences;
    }

    /**
     * Using array instead of collection provide performance improvement due to low cost of iteration and random access by index
     * @param fill
     */
    public static int[] fillPossibleSequencesFromArray(int[] fill) {
        int[] tmp = new int[fill.length];
        int last = 0;
        for (int j = 0; j < fill.length; j++) {
            if (fill[j] != 0) {
                tmp[last++] = j;
            }
        }
        int[] result = new int[last];
        System.arraycopy(tmp, 0, result, 0, last);
        return result;
    }
}
