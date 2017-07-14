package by.bsu.util;

import by.bsu.distance.LevenshteinDistance;
import by.bsu.model.Sample;
import by.bsu.model.SequencesTree;
import com.carrotsearch.hppc.LongHashSet;
import com.carrotsearch.hppc.LongSet;
import com.carrotsearch.hppc.ShortArrayList;
import info.debatty.java.stringsimilarity.QGram;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

/**
 * Just class with some useful functions
 */
public class Utils {
    private static final List<String> impossibleCharacters = new ArrayList<>();

    public static List<String> numbers = new ArrayList<>();

    static {
        for (int i = 0; i < 100_000; i++) {
            numbers.add(String.valueOf(i));
        }
    }
    static {
        impossibleCharacters.add("");
        for (int i = 1; i < 200; i++) {
            impossibleCharacters.add(impossibleCharacters.get(i-1)+"%");
        }
    }

    public static long getHashValue(int position, int l, String str) {
        long hashValue = 0;
        for (int j = 0; j < l; j++) {
            hashValue *= 4;
            hashValue += convertLetterToDigit(str.charAt(position+j));
        }
        return hashValue;
    }

    /**
     * Return hashes of l-mers in string
     */
    public static LongSet getSequenceHashesSet(String seq, int l){
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
    public static List<Long> getSequenceHashesList(String seq, int l){
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
    public static long[] getSequenceHashesArray(String seq, int l){
        long[] result = new long[seq.length() - l+1];
        int i = 0;
        long hashValue = 0;

        for (int j = 0; j < l; j++) {
            hashValue *= 4;
            hashValue += convertLetterToDigit(seq.charAt(j));
        }
        result[i++] = hashValue;
        for (int j = 1; j < seq.length() - l +1; j++) {
            hashValue -= convertLetterToDigit(seq.charAt(j-1)) << 2 * (l-1);
            hashValue <<= 2;
            hashValue += convertLetterToDigit(seq.charAt(j+l -1));
            result[i++] = hashValue;
        }
        return result;
    }

    public static void fillNodeGramsAndChunks(SequencesTree.Node node, int l){
        long hashValue = 0;
        String seq = node.key;
        if (seq.length() < l){
            return;
        }
        for (int j = 0; j < l; j++) {
            hashValue *= 4;
            hashValue += convertLetterToDigit(seq.charAt(j));
        }
        if (!node.grams.containsKey(hashValue)){
            node.grams.put(hashValue, new ShortArrayList());
        }
        node.grams.get(hashValue).add((short) 0);
        node.chunks.add(hashValue);
        for (short j = 1; j < seq.length() - l +1; j++) {
            hashValue -= convertLetterToDigit(seq.charAt(j-1)) << 2 * (l-1);
            hashValue <<= 2;
            hashValue += convertLetterToDigit(seq.charAt(j+l -1));
            if (!node.grams.containsKey(hashValue)){
                node.grams.put(hashValue, new ShortArrayList());
            }
            node.grams.get(hashValue).add(j);
            if (j % l == 0){
                node.chunks.add(hashValue);
            }
        }
    }

    /**
     * The number of hits between sequences' q-grams
     */
    public static int qHits(String s1, String s2, int q){
        long[] h1 = getSequenceHashesArray(s1, q);
        long[] h2 = getSequenceHashesArray(s2, q);
        Arrays.sort(h1);
        Arrays.sort(h2);
        return qHits(h1, h2);
    }

    public static int qHits(long[] h1, long[] h2){
       int hits = 0;
       int i = 0, j = 0;
       while (i < h1.length && j < h2.length){
           if (h1[i] == h2[j]){
               hits++;
               i++;
               j++;
           } else if (h1[i] < h2[j]){
               i++;
           } else{
               j++;
           }
       }
       return hits;
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

    public static char convertIntToLetter(int i){
        switch (i){
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
        }
        return '_';
    }


    public static String consensus(String[] sequences){
        if (sequences.length == 0){
            return "";
        }
        int l = sequences[0].length();
        int[][] count = new int[4][l];
        for (String s : sequences) {
            for (int i = 0; i < s.length(); i++) {
                count[convertLetterToDigit(s.charAt(i))][i]++;
            }
        }
        StringBuilder str = new StringBuilder();
        for (int i = 0; i < l; i++) {
            int max = 0;
            for (int j = 1; j < 4; j++) {
                if (count[j][i] > count[max][i]){
                    max = j;
                }
            }
            str.append(convertIntToLetter(max));
        }
        return str.toString();
    }

    public static String consensus(Sample sample){
        return consensus(sample.sequences);
    }

    /**
     * Append missing characters to string so they have the same size
     */
    public static String[] stringsForHamming(String[] sequences){
        int max = Arrays.stream(sequences).mapToInt(String::length).max().getAsInt();
        String[] result = new String[sequences.length];
        for (int i = 0; i < sequences.length; i++) {
            result[i] = sequences[i].length() < max? 
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
}
