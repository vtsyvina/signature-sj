package by.bsu.model;

import com.carrotsearch.hppc.IntSet;

import java.util.Map;
import java.util.Objects;

/**
 *
 * Data container to store point info for each sequence from sample
 */
public class Points {
    public String sampleName;
    /**
     * Map sequence -> coordinates where coordinates four elements array:
     * 0 -> amount of A in a certain sequence
     * 1 -> amount of C
     * 2 -> amount of G
     * 3 -> amount of T
     */
    public Map<Point, IntSet> pointSeqMap;

    public Points() {
    }

    public Points(String sampleName, Map<Point, IntSet> sequencePoints) {
        this.sampleName = sampleName;
        this.pointSeqMap = sequencePoints;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Points points = (Points) o;
        return Objects.equals(sampleName, points.sampleName) &&
                Objects.equals(pointSeqMap, points.pointSeqMap);
    }

    @Override
    public int hashCode() {
        return Objects.hash(sampleName, pointSeqMap);
    }

    @Override
    public String toString() {
        return "Points{" +
                "sampleName='" + sampleName + '\'' +
                ", pointSeqMap=" + pointSeqMap +
                '}';
    }
}
