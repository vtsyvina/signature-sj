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

    public static Map<Integer, String> readList(Path filePath) throws IOException {
        List<String> raw = Files.readAllLines(filePath);
        Map<Integer, String> result = new HashMap<>();
        int i = 0;
        StringBuilder seq = new StringBuilder();
        for (int j = 0; j < raw.size(); j++) {
            String str = raw.get(j);
            if (str.startsWith(">") && seq.length() > 0 || str.length() == 0) {
                result.put(i, seq.toString());
                i++;
                seq.setLength(0);
            } else if (!str.startsWith(">")) {
                seq.append(str);
                //workaround for for cases when file doesn't cantain last empty string
                if (j + 1 == raw.size()) {
                    result.put(i, seq.toString());
                }
            }
        }
        return result;
    }

    public static Map<Integer, String> readList(String filePath) throws IOException {
        return readList(Paths.get(filePath));
    }

    public static List<Sample> readSampleList(String filePath, boolean onePerFile) throws IOException {
        File dir = new File(filePath);
        List<Sample> result = new ArrayList<>();
        for (File file : dir.listFiles()) {
            if (onePerFile == file.isDirectory()) {
                continue;
            }
            result.add(onePerFile ? readSampleFromFile(file) : readSampleFromFolder(file));
        }
        return result;
    }

    public static Sample readSampleFromFolder(File file) {
        List<File> sortedFiles = Arrays.stream(file.listFiles())
                .sorted(Comparator.comparing(File::getName))
                .collect(Collectors.toList());
        Map<Integer, String> seq = new HashMap<>();
        sortedFiles.stream().filter( f -> !f.isHidden()).forEach(f -> {
            try {
                seq.putAll(readList(f.toPath())
                        .entrySet()
                        .stream()
                        .collect(Collectors.toMap(e -> e.getKey() + seq.size(), Map.Entry::getValue)));
            } catch (IOException e) {
                System.err.println("Error reading file for multyfile sample");
                e.printStackTrace();
            }
        });
        return new Sample(file.getName(), seq);
    }

    public static Sample readSampleFromFile(File file) throws IOException {
        return new Sample(file.getName(), readList(file.toPath()));
    }

}
