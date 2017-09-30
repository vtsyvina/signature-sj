package by.bsu.util.builders;

import by.bsu.model.SNVStructure;
import by.bsu.model.Sample;
import com.carrotsearch.hppc.IntArrayList;

public class SNVStructureBuilder {

    public static SNVStructure build(Sample sample, Sample src, double[][] srcProfile) {
        SNVStructure result = new SNVStructure();
        IntArrayList[] rows = new IntArrayList[sample.sequences[0].length()];
        result.rowMinors = new int[sample.sequences[0].length()][];

        IntArrayList[] cols = new IntArrayList[sample.sequences.length];
        result.colMinors = new int[sample.sequences.length][];
        result.rowN = new int[src.sequences[0].length()][];

        IntArrayList[] rowN = new IntArrayList[src.sequences[0].length()];
        for (int i = 0; i < src.sequences[0].length(); i++) {
            rowN[i] = new IntArrayList();
        }
        for (int i = 0; i < src.sequences.length; i++) {
            for (int j = 0; j < src.sequences[i].length(); j++) {
                if (src.sequences[i].charAt(j) == 'N'){
                    rowN[j].add(i);
                }
            }
        }

        for (int i = 0; i < sample.sequences.length; i++) {
            cols[i] = new IntArrayList();
        }
        for (int j = 0; j < sample.sequences[0].length(); j++) {
            rows[j] = new IntArrayList();
        }
        for (int i = 0; i < sample.sequences.length; i++) {
            for (int j = 0; j < sample.sequences[i].length(); j++) {
                if (sample.sequences[i].charAt(j) == '2') {
                    rows[j].add(i);
                    cols[i].add(j);
                }
            }
        }
        for (int i = 0; i < src.sequences[0].length(); i++) {
            result.rowN[i] = rowN[i].toArray();
        }
        for (int i = 0; i < sample.sequences.length; i++) {
            result.colMinors[i] = cols[i].toArray();
        }
        for (int i = 0; i < sample.sequences[0].length(); i++) {
            result.rowMinors[i] = rows[i].toArray();
        }
        result.profile = srcProfile;
        return result;
    }
}
