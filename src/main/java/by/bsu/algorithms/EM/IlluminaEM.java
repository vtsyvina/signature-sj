package by.bsu.algorithms.EM;

import by.bsu.model.IlluminaSNVSample;
import by.bsu.model.PairEndRead;

import java.util.List;

public class IlluminaEM extends AbstractEM {

    private static final double e = 0.001;

    @Override
    protected double getE() {
        return e;
    }

    public List<Double> frequencies(List<String> haplotypes, IlluminaSNVSample sample) {
        double[] frequencies = new double[haplotypes.size()];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = 1 / (double) frequencies.length;
        }
        double[][] h = new double[haplotypes.size()][sample.reads.size()];
        int f = 0;
        for (String haplotype : haplotypes) {
            int s = 0;
            for (PairEndRead read : sample.reads) {
                int misses = 0;
                for (int i = 0; i < read.l.length(); i++) {
                    if (haplotype.charAt(read.lOffset + i) != read.l.charAt(i)) {
                        misses++;
                    }
                }
                for (int i = 0; i < read.r.length(); i++) {
                    if (haplotype.charAt(read.rOffset + i) != read.r.charAt(i)) {
                        misses++;
                    }
                }
                h[f][s++] = pre[misses];
            }
            f++;
        }
        return calculateFrequencies(haplotypes.size(), sample.reads.size(), h, false);
    }
}
