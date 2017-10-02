package by.bsu.util.builders;

import by.bsu.model.IlluminaSNVSample;
import by.bsu.model.PairEndRead;
import by.bsu.model.SNVStructure;
import by.bsu.model.Sample;
import com.carrotsearch.hppc.IntArrayList;

public class SNVStructureBuilder {

    /**
     * Builds data structure for PacBio reads
     *
     * @param sample     splitted sample for all minors
     * @param src        source sample
     * @param srcProfile source profile
     * @return SNV structure with lists of minor positions
     */
    public static SNVStructure buildPacBio(Sample sample, Sample src, double[][] srcProfile) {
        SNVStructure result = new SNVStructure();

        IntArrayList[] rows = initIntArrayList(sample.sequences[0].length());
        result.rowMinors = new int[sample.sequences[0].length()][];

        IntArrayList[] cols = initIntArrayList(sample.sequences.length);
        result.colMinors = new int[sample.sequences.length][];

        result.rowN = new int[src.sequences[0].length()][];
        IntArrayList[] rowN = initIntArrayList(src.sequences[0].length());
        for (int i = 0; i < src.sequences.length; i++) {
            fillRowN(src.sequences[i], 0, i, rowN);
        }

        for (int i = 0; i < sample.sequences.length; i++) {
            fillRowsAndCols(sample.sequences[i], 0, i, cols, rows);
        }
        result.rowN = copyFromIntListToArray(rowN);
        result.colMinors = copyFromIntListToArray(cols);
        result.rowMinors = copyFromIntListToArray(rows);
        result.profile = srcProfile;
        return result;
    }

    public static SNVStructure buildIllumina(IlluminaSNVSample sample, IlluminaSNVSample src, double[][] srcProfile) {
        SNVStructure result = new SNVStructure();
        IntArrayList[] rows = initIntArrayList(sample.referenceLength);
        result.rowMinors = new int[sample.referenceLength][];

        IntArrayList[] cols = initIntArrayList(sample.reads.size());
        result.colMinors = new int[sample.reads.size()][];

        result.rowN = new int[src.referenceLength][];
        IntArrayList[] rowN = initIntArrayList(src.referenceLength);
        for (int i = 0; i < src.reads.size(); i++) {
            PairEndRead read = src.reads.get(i);
            fillRowN(read.l, read.lOffset, i, rowN);
            fillRowN(read.r, read.rOffset, i, rowN);
        }
        for (int i = 0; i < sample.reads.size(); i++) {
            PairEndRead read = sample.reads.get(i);
            fillRowsAndCols(read.l, read.lOffset, i, cols, rows);
            fillRowsAndCols(read.r, read.rOffset, i, cols, rows);
        }

        result.rowN = copyFromIntListToArray(rowN);
        result.colMinors = copyFromIntListToArray(cols);
        result.rowMinors = copyFromIntListToArray(rows);
        result.profile = srcProfile;
        return result;
    }

    private static void fillRowsAndCols(String read, int offset, int readNumber, IntArrayList[] cols, IntArrayList[] rows) {
        for (int j = 0; j < read.length(); j++) {
            if (read.charAt(j) == '2') {
                cols[readNumber].add(offset + j);
                rows[offset + j].add(readNumber);
            }
        }
    }

    private static void fillRowN(String read, int offset, int readNumber, IntArrayList[] rowN) {
        for (int j = 0; j < read.length(); j++) {
            if (read.charAt(j) == 'N') {
                rowN[readNumber].add(j + offset);
            }
        }
    }

    private static IntArrayList[] initIntArrayList(int l){
        IntArrayList[] result = new IntArrayList[l];
        for (int i = 0; i < l; i++) {
            result[i] = new IntArrayList();
        }
        return result;
    }

    private static int[][] copyFromIntListToArray(IntArrayList[] src){
        int[][] result = new int[src.length][];
        for (int i = 0; i < src.length; i++) {
            result[i] = src[i].toArray();
        }
        return result;
    }
}
