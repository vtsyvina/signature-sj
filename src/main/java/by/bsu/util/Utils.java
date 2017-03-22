package by.bsu.util;

import com.carrotsearch.hppc.LongHashSet;
import com.carrotsearch.hppc.LongSet;

/**
 * Just class with some useful functions
 */
public class Utils {

    public static LongSet getSequenceHashes(String seq, int l){
        LongSet result = new LongHashSet(seq.length() - l);
        long hashValue = 0;

        for (int j = 0; j < l; j++) {
            hashValue *= 4;
            hashValue += convertLetterToDigit(seq.charAt(j));
        }
        result.add(hashValue);
        for (int j = 1; j < seq.length() - l +1; j++) {
            hashValue -= convertLetterToDigit(seq.charAt(j-1)) << 2 * (l-1);
            hashValue <<= 2;
            hashValue += convertLetterToDigit(seq.charAt(j+l -1));
            result.add(hashValue);
        }
        return result;
    }

    public static int convertLetterToDigit(char c){
        switch (c) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
        }
        return -1;
    }

}
