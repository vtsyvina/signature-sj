package by.bsu.algorithms.EM;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public abstract class AbstractEM {

    private final double e = getE();
    protected static final double eps = 0.0001;

    protected double[] pre = new double[1000];

    {
        for (int i = 0; i < 1000; i++) {
            pre[i] = 1;
            for (int j = 0; j < 1000; j++) {
                if (j < i) {
                    pre[i] *= e / 3;
                } else {
                    pre[i] *= 1 - e;
                }
            }
        }
    }

    protected abstract double getE();


    protected double euclidDistance(double[] x, double[] y) {
        double d = 0;
        for (int i = 0; i < x.length; i++) {
            d += (x[i] - y[i]) * (x[i] - y[i]);
        }
        return Math.sqrt(d);
    }

    protected List<Double> calculateFrequencies(int haplotypesCount, int readsCount, double[][] h, boolean log) {
        double[] frequencies = new double[haplotypesCount];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = 1 / (double) frequencies.length;
        }
        double[] oldFrequencies;
        do {
            oldFrequencies = frequencies;
            frequencies = new double[haplotypesCount];
            double[] denominators = new double[readsCount];
            for (int i = 0; i < readsCount; i++) {
                double d = 0;
                for (int j = 0; j < haplotypesCount; j++) {
                    d += oldFrequencies[j] * h[j][i];
                }
                denominators[i] = d;
            }
            double[] m = new double[haplotypesCount];
            double sum = 0;
            for (int j = 0; j < haplotypesCount; j++) {

                for (int i = 0; i < readsCount; i++) {
                    if (denominators[i] == 0.0){
                        continue;
                    }
                    m[j] += oldFrequencies[j] * h[j][i] / denominators[i];
                }
                sum += m[j];
            }
            for (int j = 0; j < haplotypesCount; j++) {
                frequencies[j] = m[j] / sum;
            }
            if (log) System.out.println(Arrays.toString(frequencies));
        } while (euclidDistance(oldFrequencies, frequencies) > eps);

        return Arrays.stream(frequencies)
                .boxed()
                .collect(Collectors.toList());
    }
}
