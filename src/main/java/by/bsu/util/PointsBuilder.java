package by.bsu.util;

import by.bsu.model.Point;
import by.bsu.model.Points;
import by.bsu.model.Sample;
import com.carrotsearch.hppc.IntScatterSet;

import java.util.HashMap;
import java.util.Map;

/**
 * Util class to to build points according to sample
 */
public class PointsBuilder {

    public static Points buildPoints(Sample sample){
        Points points = new Points();
        points.sampleName = sample.name;
        points.pointSeqMap = new HashMap<>();
        for (Map.Entry<Integer, String> entry : sample.sequences.entrySet()){
            int[] coordinates = getCoordinates(entry.getValue());
            Point point = new Point(coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    coordinates[3]);
            if (!points.pointSeqMap.containsKey(point)){
                points.pointSeqMap.put(point, new IntScatterSet());
            }
            points.pointSeqMap.get(point).add(entry.getKey());
        }
        return points;
    }

    public static int[] getCoordinates(String sequence) {
        int[] coordinates = new int[4];
        for (char c : sequence.toCharArray()){
            coordinates[convertLetterToDigit(c)]++;
        }
        return coordinates;
    }

    private static int convertLetterToDigit(char c){
        switch (c) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
        }
        return -1;
    }
}
