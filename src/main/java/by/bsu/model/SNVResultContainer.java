package by.bsu.model;

import java.util.Set;


/**
 * Class contains haplotype obtained in SNV method and some additional
 * human-friendly fields to understand and analyse results
 */
public class SNVResultContainer {

    public Set<String> pacBioCluster;
    public Set<PairEndRead> illuminaCluster;
    public Clique haploClique;
    public String haplotype;
    public String clusterString;
    public Clique sourceClique;
    public double frequency;


    public SNVResultContainer(String clusterString, Set<String> pacBioCluster, Clique haplotypeClique, String haplotype) {
        this.pacBioCluster = pacBioCluster;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
    }

    public SNVResultContainer(String clusterString,  Clique haplotypeClique, String haplotype, Set<PairEndRead> illuminaCluster) {
        this.illuminaCluster = illuminaCluster;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
    }



    @Override
    public String toString() {
        return "{\n" +
                "snps=" + haploClique +
                ",\n sourse clique='" + sourceClique +
                ",\n frequency='" + frequency +
                ",\n haplotype='" + haplotype + "\n\'" +
                "}";
    }
}
