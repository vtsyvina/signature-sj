package by.bsu.algorithms;

import by.bsu.model.SNVStructure;
import by.bsu.model.Sample;
import by.bsu.util.Utils;

import java.io.IOException;

public class SNVMethod {
    static String al = "ACGT-";

    public long run(Sample sample, SNVStructure struct) throws IOException {
        double threshold = 0.5;
        int count = 0;
        for (int i = 0; i < sample.sequences.length; i++) {
            int[] hits = new int[sample.sequences.length];
            int l = struct.colMinors[i].length;
            if (l < 10) {
                continue;
            }
            for (int j = 0; j < l; j++) {
                int[] row = struct.rowMinors[struct.colMinors[i][j]];
                for (int k = 0; k < row.length; k++) {
                    hits[row[k]]++;
                }
            }

            for (int j = 0; j < hits.length; j++) {
                if (i != j) {
                    int first = i / 4;
                    int second = j / 4;
                    int allele1 = i % 4 >= Utils.getMajorAllele(struct.profile, first) ? i % 4 + 1 : i % 4;
                    int allele2 = j % 4 >= Utils.getMajorAllele(struct.profile, second) ? j % 4 + 1 : j % 4;
                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    double pij = hits[j] / (double)struct.colMinors[j].length;
                    double pi = l / (double)struct.nCount[first];
                    double increase = pij / pi;
                    if (increase > 40) {
                        System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f increase=%.2f",
                                first, second, m1, m2, l, struct.colMinors[j].length, 100 * hits[j] / (double) l, increase));
                        count++;
                    }

                }
            }
        }
        System.out.println(count);
        return 0;
    }


}
