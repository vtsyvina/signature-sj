package by.bsu.algorithms.EM;

import by.bsu.model.Sample;

import java.util.List;

public class PacBioEM extends AbstractEM {

    private static final double e = 0.01;

    public List<Double> frequencies(List<String> haplotypes, Sample sample) {
        double[][] h = new double[haplotypes.size()][sample.sequences.length];
        int f = 0;
        for (String haplotype : haplotypes) {
            int s = 0;
            for (String read : sample.sequences) {
                int misses = 0;
                for (int i = 0; i < read.length(); i++) {
                    if (haplotype.charAt(i) != read.charAt(i) && read.charAt(i) != 'N') {
                        misses++;
                    }
                }
                h[f][s++] = pre[misses];
            }
            f++;
        }
        return calculateFrequencies(haplotypes.size(), sample.sequences.length, h, false);
    }

    @Override
    protected double getE() {
        return e;
    }
}
