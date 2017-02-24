package by.bsu.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class to read input data from standard input files
 */
public class FasReader {

    public static Map<Integer, String> readList(Path filePath) throws IOException {
        List<String> raw =  Files.readAllLines(filePath);
        Map<Integer, String> result = new HashMap<>();
        for (int i = 0; i < raw.size() / 2; i++) {
            result.put(i, raw.get(i*2+1));
        }
        return result;
    }

    public static Map<Integer, String> readList(String filePath) throws IOException {
        return readList(Paths.get(filePath));
    }
}
