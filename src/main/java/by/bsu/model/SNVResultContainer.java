package by.bsu.model;

import java.util.Set;


/**
 * Class contains haplotype obtained in SNV method and some additional
 * human-friendly fields to understand and analyse results
 */
public class SNVResultContainer {

    public Set<String> cluster;
    public Clique haploClique;
    public String haplotype;
    public String clusterString;
    public Clique sourceClique;


    public SNVResultContainer(String clusterString, Set<String> cluster, Clique haplotypeClique, String haplotype) {
        this.cluster = cluster;
        this.haploClique = haplotypeClique;
        this.haplotype = haplotype;
        this.clusterString = clusterString;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SNVResultContainer container = (SNVResultContainer) o;

        if (cluster != null ? !cluster.equals(container.cluster) : container.cluster != null) return false;
        if (haploClique != null ? !haploClique.equals(container.haploClique) : container.haploClique != null)
            return false;
        if (haplotype != null ? !haplotype.equals(container.haplotype) : container.haplotype != null) return false;
        if (clusterString != null ? !clusterString.equals(container.clusterString) : container.clusterString != null)
            return false;
        return sourceClique != null ? sourceClique.equals(container.sourceClique) : container.sourceClique == null;
    }

    @Override
    public int hashCode() {
        int result = cluster != null ? cluster.hashCode() : 0;
        result = 31 * result + (haploClique != null ? haploClique.hashCode() : 0);
        result = 31 * result + (haplotype != null ? haplotype.hashCode() : 0);
        result = 31 * result + (clusterString != null ? clusterString.hashCode() : 0);
        result = 31 * result + (sourceClique != null ? sourceClique.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "{\n" +
                "snvs=" + haploClique +
                ",\n haplotype='" + haplotype + "\n\'" +
                "}";
    }
}
