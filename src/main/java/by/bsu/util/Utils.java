package by.bsu.util;

import by.bsu.model.Sample;
import com.carrotsearch.hppc.LongHashSet;
import com.carrotsearch.hppc.LongSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Just class with some useful functions
 */
public class Utils {
    private static final List<String> impossibleCharacters = new ArrayList<>();
    private static final char impossibleChar = '%';
    public static final String DEFAULT_ALPHABET = "ACGT";

    public static List<String> numbers = new ArrayList<>();

    static {
        for (int i = 0; i < 100_000; i++) {
            numbers.add(String.valueOf(i));
        }
    }

    static {
        impossibleCharacters.add("");
        for (int i = 1; i < 3000; i++) {
            impossibleCharacters.add(impossibleCharacters.get(i - 1) + impossibleChar);
        }
    }

    public static long getHashValue(int position, int l, String str) {
        return getHashValue(position, l, str, DEFAULT_ALPHABET);
    }

    public static long getHashValue(int position, int l, String str, String alphabet) {
        long hashValue = 0;
        for (int j = 0; j < l; j++) {
            hashValue *= alphabet.length();
            hashValue += convertLetterToDigit(str.charAt(position + j), alphabet);
        }
        return hashValue;
    }

    /**
     * Return hashes of l-mers in string
     */
    public static LongSet getSequenceHashesSet(String seq, int l) {
        LongSet result = new LongHashSet(seq.length() - l + 1);
        long[] sequenceHashesArray = getSequenceHashesArray(seq, l);
        for (long l1 : sequenceHashesArray) {
            result.add(l1);
        }
        return result;
    }

    /**
     * Return hashes of l-mers in string
     */
    public static List<Long> getSequenceHashesList(String seq, int l) {
        List<Long> result = new ArrayList<>(seq.length() - l + 1);
        long[] sequenceHashesArray = getSequenceHashesArray(seq, l);
        for (long l1 : sequenceHashesArray) {
            result.add(l1);
        }
        return result;
    }

    /**
     * Return hashes of l-mers in string
     */
    public static long[] getSequenceHashesArray(String seq, int l) {
        long[] result = new long[seq.length() - l + 1];
        int i = 0;
        long hashValue = 0;

        for (int j = 0; j < l; j++) {
            hashValue *= 4;
            hashValue += convertLetterToDigit(seq.charAt(j));
        }
        result[i++] = hashValue;
        for (int j = 1; j < seq.length() - l + 1; j++) {
            hashValue -= convertLetterToDigit(seq.charAt(j - 1)) << 2 * (l - 1);
            hashValue <<= 2;
            hashValue += convertLetterToDigit(seq.charAt(j + l - 1));
            result[i++] = hashValue;
        }
        return result;
    }


    public static int convertLetterToDigit(char c) {
        return convertLetterToDigit(c, DEFAULT_ALPHABET);
    }

    public static int convertLetterToDigit(char c, String alphabet) {
        return alphabet.indexOf(c);
    }

    public static char convertIntToLetter(int i) {
        switch (i) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            case 4:
                return '-';
        }
        return '_';
    }

    public static String consensus(String[] sequences){
        return consensus(sequences, DEFAULT_ALPHABET);
    }

    public static String consensus(String[] sequences, String alphabet) {
        if (sequences.length == 0) {
            return "";
        }
        int l = sequences[0].length();
        int[][] count = new int[alphabet.length()][l];
        for (String s : sequences) {
            for (int i = 0; i < s.length(); i++) {
                count[convertLetterToDigit(s.charAt(i), alphabet)][i]++;
            }
        }
        StringBuilder str = new StringBuilder();
        int alphabetLength = alphabet.indexOf('N') == -1 ? alphabet.length() : alphabet.length() - 1;
        for (int i = 0; i < l; i++) {
            int max = 0;
            for (int j = 1; j < alphabetLength; j++) {
                if (count[j][i] > count[max][i]) {
                    max = j;
                }
            }
            str.append(convertIntToLetter(max));
        }
        return str.toString();
    }

    public static String consensus(Sample sample) {
        return consensus(sample.sequences);
    }

    /**
     * Return consensus for given profile. Put $ if all profile values equals 0
     * @param profile sample profile
     * @param alphabet alphabet for profile
     * @return consensus
     */
    public static String consensus(double[][] profile, String alphabet){
        char[] majors = new char[profile[0].length];
        for (int j = 0; j < profile[0].length; j++) {
            int majorAllele = Utils.getMajorAllele(profile, j);
            majors[j] = majorAllele == -1 ? '$' :  alphabet.charAt(majorAllele);
        }
        return new String(majors);
    }

    public static double[][] profile(Sample sample){
        return profile(sample, DEFAULT_ALPHABET);
    }
    public static double[][] profile(Sample sample, String alphabet) {
        String[] sequences = sample.forHamming;
        if (sequences.length == 0) {
            return new double[0][alphabet.length()];
        }
        int l = sequences[0].length();
        int[][] count = new int[alphabet.length()][l];
        for (String s : sequences) {
            for (int i = 0; i < s.length(); i++) {
                int d = convertLetterToDigit(s.charAt(i), alphabet);
                if (d != -1) {
                    count[d][i]++;
                }
            }
        }
        double[][] result = new double[alphabet.length()][l];
        for (int i = 0; i < alphabet.length(); i++) {
            for (int j = 0; j < l; j++) {
                result[i][j] = count[i][j] / (double) sequences.length;
            }
        }
        return result;
    }


    /**
     * Append missing characters to string so they have the same size
     */
    public static String[] stringsForHamming(String[] sequences) {
        int max = Arrays.stream(sequences).mapToInt(String::length).max().getAsInt();
        String[] result = new String[sequences.length];
        for (int i = 0; i < sequences.length; i++) {
            result[i] = sequences[i].length() < max ?
                    sequences[i] + impossibleCharacters.get(max - sequences[i].length()) :
                    sequences[i];
        }
        return result;
    }

    public static void expandNumbers(int size) {
        if (size > numbers.size()) {
            for (int i = numbers.size(); i < size; i++) {
                numbers.add(String.valueOf(i));
            }
        }
    }

    public static int getMajorAllele(double[][] profile, int i) {
        int major = 0;
        for (int j = 1; j < profile.length; j++) {
            if (profile[j][i] > profile[major][i]){
                major = j;
            }
        }
        //in case there is no reads covering this position
        return profile[major][i] < 0.001 ? -1 : major;
    }

}
