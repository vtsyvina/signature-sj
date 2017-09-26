package by.bsu.model;

public class SNVSample extends Sample {
    public int[] offset;
    public int[] readLength;

    public SNVSample(String name, String[] sequences) {
        super(name, sequences);
    }

    public SNVSample(String name, String[] sequences, int[] offset, int[] readLength) {
        super(name, sequences);
        this.offset = offset;
        this.readLength = readLength;
    }
}
