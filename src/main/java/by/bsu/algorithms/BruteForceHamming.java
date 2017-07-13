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
        for (Map.Entry<Integer, String> entry1 : sample1.sequences.entrySet()){
            iter++;
            if (iter % 400 == 0){
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r" + iter);
            }
            for (Map.Entry<Integer, String> entry2 : sample2.sequences.entrySet()){
                int d = distance.apply(entry1.getValue(), entry2.getValue());
                if (d <= k){
                    result++;
                    str.append(numbers.get(entry1.getKey())).append(" ").append(numbers.get(entry2.getKey())).append("\n");
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
        int i = 0;
        Path path = Start.getOutputFilename(sample, "brute-hamming");
        StringBuilder str = new StringBuilder();
        for (Map.Entry<Integer, String> entry1 : sample.sequences.entrySet()){
            i++;
            if (i % 400 == 0){
                Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
                str = new StringBuilder();
                System.out.print("\r"+i);
            }
            for (Map.Entry<Integer, String> entry2 : sample.sequences.entrySet()){
                if (entry1.getKey() <= entry2.getKey()){
                    continue;
                }
                int d = distance.apply(entry1.getValue(), entry2.getValue());
                if (d <= k){
                    result++;
                    str.append(numbers.get(entry1.getKey())).append(" ").append(numbers.get(entry2.getKey())).append("\n");
                }
            }
        }
        Files.write(path, str.toString().getBytes(), StandardOpenOption.APPEND);
        System.out.println();
        System.out.println("length = " +result);
        return result;
    }
}
