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
        for (LongSet positionHashes : dict1.wholeSampleFixedPositionHashesList) {
            for (LongCursor hash : positionHashes) {
                if (dict2.allHashesSet.contains(hash.value)) {
                    kMerCoincidences++;
                    break;
                }
            }
        }
        return kMerCoincidences;
    }
}
