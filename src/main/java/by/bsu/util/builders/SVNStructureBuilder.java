package by.bsu.util.builders;

import by.bsu.model.SVNStructure;
import by.bsu.model.Sample;
import com.carrotsearch.hppc.IntArrayList;

import java.util.Arrays;

public class SVNStructureBuilder {

    public static SVNStructure build(Sample rotated, Sample src, double[][] srcProfile) {
        SVNStructure result = new SVNStructure();
        IntArrayList[] rows = new IntArrayList[rotated.sequences[0].length()];
        result.rowMinors = new int[rotated.sequences[0].length()][];
        IntArrayList[] cols = new IntArrayList[rotated.sequences.length];
        result.colMinors = new int[rotated.sequences.length][];
        result.nCount = new int[src.sequences.length];
        for (int i = 0; i < src.sequences.length; i++) {
            result.nCount[i] = src.sequences.length;
        }
        for (String sequence : src.sequences) {
            for (int i = 0; i < sequence.length(); i++) {
                if (sequence.charAt(i) == 'N'){
                    result.nCount[i]--;
                }
            }
        }

        for (int i = 0; i < rotated.sequences.length; i++) {
            cols[i] = new IntArrayList();
        }
        for (int j = 0; j < rotated.sequences[0].length(); j++) {
            rows[j] = new IntArrayList();
        }
        for (int i = 0; i < rotated.sequences.length; i++) {
            for (int j = 0; j < rotated.sequences[i].length(); j++) {
                if (rotated.sequences[i].charAt(j) == '2') {
                    rows[j].add(i);
                    cols[i].add(j);
                }
            }
        }
        for (int i = 0; i < rotated.sequences.length; i++) {
            result.colMinors[i] = cols[i].toArray();
        }
        for (int i = 0; i < rotated.sequences[0].length(); i++) {
            result.rowMinors[i] = rows[i].toArray();
        }
        result.profile = srcProfile;
        return result;
    }
}
