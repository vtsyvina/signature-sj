package by.bsu.util;

import by.bsu.model.SVNSample;
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

    public static SVNSample readOneLined(File file) throws IOException{
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

        return new SVNSample(file.getName(),
                Utils.stringsForHamming(tmp.toArray(new String[0]), 'N'),
                offset.stream().mapToInt(o-> o).toArray(), readLength.stream().mapToInt(o-> o).toArray());
    }

    public static Sample rowsColsRotate(Sample sample, double[][] profile, String alphabet){
        int minorsCount = alphabet.length()-1;
        String[] sequences = new String[sample.sequences[0].length()*minorsCount];
        StringBuilder str = new StringBuilder();
        for (int i = 0; i < sample.sequences[0].length(); i++) {
            int major = Utils.getMajorAllele(profile, i);//current major
            int minor = 0; //current minor
            for (int j = 0; j < minorsCount; j++, minor++) { //change current minor + skip major, so we will have 4 rows for 4 different minors from 'ACGT-'
                if (minor == major){
                    minor++;
                }
                str.setLength(0);
                for (int k = 0; k < sample.sequences.length; k++) {
                    int allele = alphabet.indexOf(sample.sequences[k].charAt(i));
                    if (allele == minor){
                        str.append("2");
                    } else {
                        str.append("1");
                    }
                }
                sequences[i*minorsCount+j] = str.toString();
            }
        }
        return new Sample(sample.name+"_rotated", sequences);
    }

}
