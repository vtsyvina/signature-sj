package by.bsu.util;

import by.bsu.model.SNVSample;
import by.bsu.model.Sample;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class to read input data from standard input files
 */
public class FasReader {

    public static String[] readList(Path filePath) throws IOException {
        return readList(filePath, false);
    }

    public static String[] readList(Path filePath, boolean reversed) throws IOException{
        List<String> raw = Files.readAllLines(filePath);
        List<String> result = new ArrayList<>();
        int i = 0;
        StringBuilder seq = new StringBuilder();
        for (int j = 0; j < raw.size(); j++) {
            String str = raw.get(j);
            if (str.startsWith(">") && seq.length() > 0 || str.length() == 0) {
                result.add((reversed ? seq.reverse() : seq).toString());
                i++;
                seq.setLength(0);
            } else if (!str.startsWith(">")) {
                seq.append(str);
                //workaround for for cases when file doesn't contain last empty string
                if (j + 1 == raw.size()) {
                    result.add((reversed ? seq.reverse() : seq).toString());
                }
            }
        }
        return result.stream().toArray(String[]::new);
    }

    public static String[] readList(String filePath) throws IOException {
        return readList(Paths.get(filePath));
    }

    public static List<Sample> readSampleList(String filePath, boolean onePerFile) throws IOException{
        return readSampleList(filePath, onePerFile, false);
    }

    public static List<Sample> readSampleList(String filePath, boolean onePerFile, boolean reversed) throws IOException {
        File dir = new File(filePath);
        List<Sample> result = new ArrayList<>();
        for (File file : dir.listFiles()) {
            if (onePerFile == file.isDirectory()) {
                continue;
            }
            result.add(onePerFile ? readSampleFromFile(file) : readSampleFromFolder(file, reversed));
        }
        return result;
    }

    public static Sample readSampleFromFolder(File file) {
        return readSampleFromFolder(file, false);
    }

    public static Sample readSampleFromFolder(File file, boolean reversed) {
        List<File> sortedFiles = Arrays.stream(file.listFiles())
                .sorted(Comparator.comparing(File::getName))
                .collect(Collectors.toList());
        List<String> seq = new ArrayList<>();
        sortedFiles.stream().filter( f -> !f.isHidden()).forEach(f -> {
            try {
                seq.addAll(Arrays.asList(readList(f.toPath(), reversed)));
            } catch (IOException e) {
                System.err.println("Error reading file for multyfile sample");
                e.printStackTrace();
            }
        });
        
        return new Sample(file.getName(), seq.stream().toArray(String[]::new));
    }

    public static Sample readSampleFromFile(File file) throws IOException {
        return new Sample(file.getName(), readList(file.toPath()));
    }

    public static SNVSample readOneLined(File file) throws IOException{
        List<String> raw = Files.readAllLines(file.toPath());
        List<String> tmp = new ArrayList<>();
        List<String> N = new ArrayList<>();
        List<Integer> offset = new ArrayList<>();
        List<Integer> readLength = new ArrayList<>();
        N.add("");
        for (int i = 1; i < 3000; i++) {
            N.add(N.get(i - 1) + 'N');
        }
        int i = 0;
        for (String s : raw) {
            String[] split = s.split("\t");
            offset.add(Integer.valueOf(split[3])-1);
            readLength.add(split[5].length());
            tmp.add(N.get(offset.get(i))+split[5]);
            i++;
        }

        return new SNVSample(file.getName(),
                Utils.stringsForHamming(tmp.toArray(new String[0]), 'N'),
                offset.stream().mapToInt(o-> o).toArray(), readLength.stream().mapToInt(o-> o).toArray());
    }

    /**
     * Splits each column of the given sample into 4 columns for each minor. For example:
     * Given alphabet as "ACGT-"
     * for i-th column 'G' is major. then i-th column will split into 4 columns where first column will place '2' for each 'A' and '1' for any other allele,
     * for second column will place '2' for each 'C' and '1' otherwise, for third '1' for each 'T' and for forth '1' for each '-'
     * @param sample
     * @param profile
     * @param alphabet
     * @return
     */
    public static Sample splitColumns(Sample sample, double[][] profile, String alphabet){
        int minorsCount = alphabet.length()-1;
        String[] sequences = new String[sample.sequences.length];
        StringBuilder str = new StringBuilder();
        for (int i = 0; i < sample.sequences.length; i++) {
            str.setLength(0);
            for (int j = 0; j < sample.sequences[i].length(); j++) {
                int major = Utils.getMajorAllele(profile, j);
                int minor = 0;
                for (int k = 0; k < minorsCount; k++, minor++) {
                    if (minor == major){
                        minor++;
                    }
                    int allele = alphabet.indexOf(sample.sequences[i].charAt(j));
                    str.append(allele == minor ? "2" : "1");
                }
            }
            sequences[i] = str.toString();
        }

        return new Sample(sample.name+"_splitted", sequences);
    }

}
