package by.bsu.util;

import by.bsu.model.Sample;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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

}
