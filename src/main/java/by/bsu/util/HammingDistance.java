package by.bsu.util;

import org.apache.commons.text.beta.similarity.EditDistance;

/**
 * Hamming distance with threshold
 */
public class HammingDistance {

    public static Integer apply(CharSequence left, CharSequence right, int threshold) {
        if (left == null || right == null) {
            throw new IllegalArgumentException("Strings must not be null");
        }

        if (left.length() != right.length()) {
            throw new IllegalArgumentException("Strings must have the same length");
        }

        int distance = 0;

        for (int i = 0; i < left.length(); i++) {
            if (left.charAt(i) != right.charAt(i)) {
                distance++;
                if (distance > threshold){
                    return -1;
                }
            }
        }
        return distance;
    }
}
