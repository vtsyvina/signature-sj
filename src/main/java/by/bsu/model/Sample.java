package by.bsu.model;

import by.bsu.util.Utils;

import java.util.HashMap;
import java.util.Map;

/**
 * Data container to store sample data
 */
public class Sample {
    public String name;
    /**
     * Contains list of all sequences with digits instead of letters A -> 0, C -> 1, G -> 2, T -> 3
     */
    public Map<Integer, String> sequences;

    /**
     * In case we want faster performance avoiding String.charAt method
     */
    public Map<Integer, char[]> sequencesChars;

    public Map<Integer, String> forHamming;
    
    public Map<Integer, Map<String, Integer>> profiles;

    public String consensus;

    public Sample() {
    }

    public Sample(String name, Map<Integer, String> sequences) {
        this.name = name;
        this.sequences = sequences;
        this.forHamming = Utils.stringsForHamming(sequences);
        this.sequencesChars = new HashMap<>();
        sequences.entrySet().stream().forEach( e-> {
            int l = e.getValue().length();
            char[] tmp = new char[l];
            e.getValue().getChars(0, l, tmp, 0);
            sequencesChars.put(e.getKey(), tmp);
        });
    }
    
    public Sample(String name, Map<Integer, String> sequences, int l){
        this(name, sequences);
        this.profiles = Utils.getProfiles(sequences, l);
    }

    public Sample(String name, Map<Integer, String> sequences, String consensus) {
        this(name, sequences);
        this.forHamming = Utils.stringsForHamming(sequences);
    }
}
