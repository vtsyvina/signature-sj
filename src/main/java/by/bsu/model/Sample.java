package by.bsu.model;

import by.bsu.util.Utils;

/**
 * Data container to store sample data
 */
public class Sample {
    public String name;
    /**
     * Contains list of all sequences with digits instead of letters A -> 0, C -> 1, G -> 2, T -> 3
     */
    public String[] sequences;

    /**
     * In case we want faster performance avoiding String.charAt method
     */
    public char[][] sequencesChars;

    public String[] forHamming;

    public String consensus;

    public Sample() {
    }

    public Sample(String name, String[] sequences) {
        this.name = name;
        this.sequences = sequences;
        this.forHamming = Utils.stringsForHamming(sequences);
        this.sequencesChars = new char[sequences.length][];
        for (int i = 0; i < sequences.length; i++) {
            int l = sequences[i].length();
            char[] tmp = new char[l];
            sequences[i].getChars(0, l, tmp, 0);
            sequencesChars[i] = tmp;
        }
    }
    
    public Sample(String name, String[] sequences, int l){
        this(name, sequences);
    }
}
