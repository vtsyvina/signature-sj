package by.bsu.algorithms;

import by.bsu.model.SNVStructure;
import by.bsu.model.Sample;
import by.bsu.util.Utils;

import java.io.IOException;

public class SNVMethod {
    static String al = "ACGT-";

    public long run(Sample sample, SNVStructure struct, Sample src) throws IOException {
        double threshold = 0.5;
        int count = 0;
        int countP = 0;
        for (int i = 0; i < sample.sequences[0].length(); i++) {
            int[] hits = new int[sample.sequences[0].length()];
            int l = struct.rowMinors[i].length;
            if (l < 5) {
                continue;
            }
            for (int j = 0; j < l; j++) {
                int[] column = struct.colMinors[struct.rowMinors[i][j]];
                for (int aColumn : column) {
                    hits[aColumn]++;
                }
            }

            for (int j = 0; j < hits.length; j++) {
                if (i != j && hits[j] >= 10) {
                    int first = i / 4;
                    int second = j / 4;
                    int allele1 = i % 4 >= Utils.getMajorAllele(struct.profile, first) ? i % 4 + 1 : i % 4;
                    int allele2 = j % 4 >= Utils.getMajorAllele(struct.profile, second) ? j % 4 + 1 : j % 4;
                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);

                    int o22 = hits[j];
                    int o21 = struct.rowMinors[i].length-hits[j];
                    int o12 = struct.rowMinors[j].length-hits[j];

                    for (int k = 0; k < struct.rowMinors[i].length; k++) {
                        if (src.sequences[struct.rowMinors[i][k]].charAt(second) == 'N'){
                            o21--;
                        }
                    }
                    for (int k = 0; k < struct.rowMinors[j].length; k++) {
                        if (src.sequences[struct.rowMinors[j][k]].charAt(first) == 'N'){
                            o12--;
                        }
                    }

                    int o11 = sample.sequences.length - o12;
                    int minNColumn = struct.rowN[first].length < struct.rowN[second].length ? first : second;
                    int maxNColumn = struct.rowN[first].length > struct.rowN[second].length ? first : second;
                    for (int k = 0; k < struct.rowN[minNColumn].length; k++) {
                        if (src.sequences[struct.rowN[minNColumn][k]].charAt( maxNColumn) == 'N'){
                            o11--;
                        }
                    }
                    //double p = struct.rowMinors[i].length/(double)(sample.sequences.length - struct.rowN[first].length);
                    double p = (o12*o21)/((double)o11*struct.rowN[first].length);
                    if (p < 0.000_000_001){
                        countP++;
                    }
                    double pij = hits[j] / (double)struct.rowMinors[j].length;
                    double pi = l / (double)struct.rowN[first].length;
                    double increase = pij / pi;
                    if (p < 1E-12) {
                        System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.2f zero",
                                first, second, m1, m2, l, struct.rowMinors[j].length, 100 * hits[j] / (double) l, p));
                        count++;
                    }else{
                        double pvalue = 1;
                        for (int k = 0; k < o22; k++) {
                            pvalue -= Utils.poissonBinomialAppr(k, p, o11+o12+o21+o22);
                        }
                        if (pvalue < 0.01){
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.2f",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, 100 * hits[j] / (double) l, pvalue));
                        }
                        count++;
                    }

                }
            }
        }
        System.out.println(count);
        System.out.println(countP);
        return 0;
    }


}
