package by.bsu.distance;


public class HammingDistance {

    /**
     * Find the Hamming Distance between two strings with the same
     * length.
     * <p>
     * <p>The distance starts with zero, and for each occurrence of a
     * different character in either String, it increments the distance
     * by 1, and finally return its value.</p>
     * <p>
     * <p>Since the Hamming Distance can only be calculated between strings of equal length, input of different lengths
     * will throw IllegalArgumentException</p>
     * <p>
     * <pre>
     * distance.apply("", "")               = 0
     * distance.apply("pappa", "pappa")     = 0
     * distance.apply("1011101", "1011111") = 1
     * distance.apply("ATCG", "ACCC")       = 2
     * distance.apply("karolin", "kerstin"  = 3
     * </pre>
     *
     * @param left  the first CharSequence, must not be null
     * @param right the second CharSequence, must not be null
     * @return distance
     * @throws IllegalArgumentException if either input is {@code null} or
     *                                  if they do not have the same length
     */
    public int apply(final CharSequence left, final CharSequence right) {
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
            }
        }

        return distance;
    }

    /**
     * Find the Hamming Distance between two char arrays with the same
     * length.
     *
     * @param left  first string
     * @param right second string
     * @param cl    arrays with chunk hashes of first string ( see {@link by.bsu.model.KMerDict} or {@link by.bsu.model.KMerDictChunks}
     * @param cr    arrays with chunk hashes of second string
     * @param l     chunks length
     * @param k     threshold
     * @return hamming distance between two strings if it is less or equal to k, returns -1 otherwise
     */
    public int apply(char[] left, char[] right, long[] cl, long[] cr, int l, int k) {
        int distance = 0;
        for (int i = 0; i < cl.length; i++) {
            if (cl[i] != cr[i]) {
                int il = i*l;
                for (int j = 0; j < l; j++) {
                    if (left[il + j] != right[il + j]) {
                        distance++;
                        if (distance > k) {
                            return -1;
                        }
                    }
                }
            }
        }
        if (left.length % l != 0) {
            int border = left.length - left.length % l;
            for (int i = left.length - 1; i >= border; i--) {
                if (left[i] != right[i]) {
                    distance++;
                }
            }
        }
        return distance > k ? -1 : distance;
    }
}
