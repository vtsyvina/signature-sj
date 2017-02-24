package by.bsu.model;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Data container to store sample data
 */
public class Sample {
    public String name;
    /**
     * Contains list of all sequences with digits instead of letters A -> 0, C -> 1, G -> 2, T -> 3
     */
    public Map<Integer, String> sequences;

    public Sample() {
    }

    public Sample(String name, Map<Integer, String> sequences) {
        this.name = name;
        this.sequences = sequences;
    }
}
