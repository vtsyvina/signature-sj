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
        for (int i = 0; i < sample.sequences[0].length(); i++) {
            int[] hits = new int[sample.sequences[0].length()];
            int l = struct.rowMinors[i].length;
            if (l < 5) {
                continue;
            }
            for (int j = 0; j < l; j++) {
                int[] column = struct.colMinors[struct.rowMinors[i][j]];
                for (int k = 0; k < column.length; k++) {
                    hits[column[k]]++;
                }
            }

            for (int j = 0; j < hits.length; j++) {
                if (i != j && hits[j] > 0) {
                    int first = i / 4;
                    int second = j / 4;
                    int allele1 = i % 4 >= Utils.getMajorAllele(struct.profile, first) ? i % 4 + 1 : i % 4;
                    int allele2 = j % 4 >= Utils.getMajorAllele(struct.profile, second) ? j % 4 + 1 : j % 4;
                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    double pij = hits[j] / (double)struct.rowMinors[j].length;
                    double pi = l / (double)struct.readsCount[first];
                    double increase = pij / pi;
                    if (increase > 40) {
                        System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f increase=%.2f",
                                first, second, m1, m2, l, struct.rowMinors[j].length, 100 * hits[j] / (double) l, increase));
                        count++;
                    }

                }
            }
        }
        System.out.println(count);
        return 0;
    }


}
