package by.bsu.algorithms;

import by.bsu.model.IlluminaSNVSample;
import by.bsu.model.PairEndRead;
import by.bsu.model.SNVStructure;
import by.bsu.util.Utils;

public class SNVIlluminaMethod {
    private static String al = "ACGT-N";
    private static int minorCount = al.length() - 2;

    public void run(IlluminaSNVSample splittedSample, SNVStructure struct, IlluminaSNVSample src) {
        int count = 0;
        String consensus = Utils.consensus(struct.profile, al);
        for (int i = 0; i < splittedSample.referenceLength; i++) {
            int l = struct.rowMinors[i].length;
            if (l < 10) {
                continue;
            }
            int[] hits = getHits(struct, struct.rowMinors[i], splittedSample.referenceLength);

            for (int j = 0; j < hits.length; j++) {
                //skip small amount of hits
                if (i != j && hits[j] >= 10) {
                    //get unsplitted columns, minors, o_kl
                    int first = i / minorCount;
                    int second = j / minorCount;
                    int allele1 = i % minorCount >= Utils.getMajorAllele(struct.profile, first) ? i % minorCount + 1 : i % minorCount;
                    int allele2 = j % minorCount >= Utils.getMajorAllele(struct.profile, second) ? j % minorCount + 1 : j % minorCount;
                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    /*
                     * false 1 means that in actual sample it has another minor or N in given position
                     */
                    int o22 = hits[j];
                    int o21 = struct.rowMinors[i].length; //all 2*
                    int o12 = struct.rowMinors[j].length; //all *2
                    // subtract 2N and false 21 from o21
                    for (int k = 0; k < struct.rowMinors[i].length; k++) {
                        PairEndRead read = src.reads.get(struct.rowMinors[i][k]);
                        //find in what part of read is second
                        if (read.lOffset <= second && read.lOffset + read.l.length() >= second) {
                            if (read.l.charAt(second - read.lOffset) != consensus.charAt(second)) {
                                o21--;
                            }
                        } else if (read.rOffset <= second && read.rOffset + read.r.length() >= second) {
                            if (read.r.charAt(second - read.rOffset) != consensus.charAt(second)) {
                                o21--;
                            }
                        }

                    }
                    //subtract N2 and false 12 from o12
                    for (int k = 0; k < struct.rowMinors[j].length; k++) {
                        PairEndRead read = src.reads.get(struct.rowMinors[j][k]);
                        //find in what part of read is first
                        if (read.lOffset <= first && read.lOffset + read.l.length() >= first) {
                            if (read.l.charAt(first - read.lOffset) != consensus.charAt(first)) {
                                o12--;
                            }
                        } else if (read.rOffset <= first && read.rOffset + read.r.length() >= first) {
                            if (read.r.charAt(first - read.rOffset) != consensus.charAt(first)) {
                                o12--;
                            }
                        }
                    }

                    int o11 = struct.majorsInRow[first]; //all 11(and some false positive 11),12,1N
                    //subtract 1N from 011
                    for (int k = 0; k < struct.rowN[second].length; k++) {
                        PairEndRead read = src.reads.get(struct.rowN[second][k]);
                        o11 = checkTrueMajor(consensus, first, o11, read);
                    }
                    //subtract 12 from 011
                    for (int k = 0; k < minorCount; k++) {
                        int[] rowMinor = struct.rowMinors[second * minorCount + k];
                        for (int m = 0; m < rowMinor.length; m++) {
                            PairEndRead read = src.reads.get(rowMinor[m]);
                            o11 = checkTrueMajor(consensus, first, o11, read);
                        }
                    }
                    // TODO what to do actually?
                    if (o11 == 0) {
                        continue;
                    }
                    //amount of common reads for i and j column
                    int reads = o11 + o12 + o21 + o22;
                    //start calculate p-value, starting with p
                    //double p = struct.rowMinors[i].length/(double)(sample.reads.length - struct.rowN[first].length);
                    double p = (o12 * o21) / ((double) o11 * reads);
                    if (p < 1E-12) {
                        System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e zero",
                                first, second, m1, m2, l, struct.rowMinors[j].length, 100 * hits[j] / (double) l, p));
                        count++;
                    } else {

                        double pvalue = Utils.binomialPvalue(o22, p, reads);
                        if (pvalue < 0.01 / (reads * (reads - 1) / 2)) {
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, 100 * hits[j] / (double) l, pvalue));
                            count++;
                        }
                    }

                }
            }
        }
        System.out.println();
    }


    /**
     * Check that given read has major allele on given position
     *
     * @param consensus Whole sample consensus
     * @param position  Position to check
     * @param o11       current o11 value to correct
     * @param read      given read
     * @return o11 if read doesn't have major on given position, o11-1 otherwise
     */
    private int checkTrueMajor(String consensus, int position, int o11, PairEndRead read) {
        if (read.lOffset <= position && read.lOffset + read.l.length() >= position && read.l.charAt(position - read.lOffset) == consensus.charAt(position)) {
            o11--; //subtract only if read has i column and major on it
        } else if (read.rOffset <= position && read.rOffset + read.r.length() >= position && read.r.charAt(position - read.rOffset) == consensus.charAt(position)) {
            o11--;
        }
        return o11;
    }

    /**
     * Calculates hits for given row of minors.
     * Gets all minors in column for a certain position for minorRow and increments appropriate hits position
     *
     * @param struct          SNV data structure
     * @param rowMinor        array for all minors for some position
     * @param referenceLength Reference length
     * @return array with all minor hits for given row
     */
    private int[] getHits(SNVStructure struct, int[] rowMinor, int referenceLength) {
        int[] hits = new int[referenceLength];

        for (int j = 0; j < rowMinor.length; j++) {
            int[] column = struct.colMinors[rowMinor[j]];
            for (int aColumn : column) {
                hits[aColumn]++;
            }
        }
        return hits;
    }
}
