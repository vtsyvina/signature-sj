package by.bsu.model;

public class SVNSample extends Sample {
    public int[] offset;
    public int[] readLength;

    public SVNSample(String name, String[] sequences) {
        super(name, sequences);
    }

    public SVNSample(String name, String[] sequences, int[] offset, int[] readLength) {
        super(name, sequences);
        this.offset = offset;
        this.readLength = readLength;
    }
}
