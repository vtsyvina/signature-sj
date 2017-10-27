package by.bsu.model;

public class SNVStructure {
    public double[][] profile;
    public int[][] rowMinors;
    public int[][] rowN;
    public int[][] colMinors;
    public int[] majorsInRow;
    //for Illumina only
    public int[][] readsAtPosition;
}
