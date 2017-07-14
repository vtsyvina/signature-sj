package by.bsu.algorithms;

import static by.bsu.util.Utils.numbers;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Map;

import by.bsu.model.Sample;
import by.bsu.start.Start;
import by.bsu.distance.HammingDistance;

/**
 * Brute force algorithm that compares all possible sequence pairs with Hamming distance
 */
public class BruteForceHamming {

    public static long run(Sample sample1, Sample sample2, int k) throws IOException {
        System.out.println("Start Brute Hamming force method for "+sample1.name+" "+sample2.name+" k="+k);
        long result = 0;
        HammingDistance distance = new HammingDistance();
        StringBuilder str = new StringBuilder();
        int iter = 0;
        Path path = Start.getOutputFilename(sample1, sample2, "brute-hamming");
        for (int i = 0; i < sample1.sequences.length; i++) {
            if (i % 400 == 0){
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r" + iter);
            }
            for (int j = 0; j < sample2.sequences.length; j++) {
                int d = distance.apply(sample1.sequences[i], sample2.sequences[j]);
                if (d <= k){
                    result++;
                    str.append(numbers.get(i)).append(" ").append(numbers.get(j)).append("\n");
                }
            }
        }
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println("length = " +result);
        return result;
    }

    public static long run(Sample sample, int k) throws IOException {
        System.out.println("Start Brute Hamming force method for "+sample.name+" k="+k);
        long result = 0;
        HammingDistance distance = new HammingDistance();
        Path path = Start.getOutputFilename(sample, "brute-hamming");
        StringBuilder str = new StringBuilder();
        for (int i = 0; i < sample.sequences.length; i++) {
            if (i % 400 == 0){
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r"+i);
            }
            for (int j = i+1; j < sample.sequences.length; j++) {
                int d = distance.apply(sample.sequences[i], sample.sequences[j]);
                if (d <= k){
                    result++;
                    str.append(numbers.get(i)).append(" ").append(numbers.get(j)).append("\n");
                }
            }
        }
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println();
        System.out.println("length = " +result);
        return result;
    }
}
