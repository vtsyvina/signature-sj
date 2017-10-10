package by.bsu.model;

import by.bsu.algorithms.SNVPacBioMethod;
import by.bsu.util.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class Clique {
    public String minors;
    public List<Integer> snps;
    public Set<Integer> splittedSnps;

    public Clique(Set<Integer> splittedSnps, SNVStructure structure) {
        this.snps = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        splittedSnps.stream().sorted().forEach(i -> {
            int pos = i / SNVPacBioMethod.minorCount;
            int allele1 = i % SNVPacBioMethod.minorCount >= Utils.getMajorAllele(structure.profile, pos) ? i % SNVPacBioMethod.minorCount + 1 : i % SNVPacBioMethod.minorCount;
            char m1 = SNVPacBioMethod.al.charAt(allele1);
            this.snps.add(pos);
            str.append(m1);
        });
        this.minors = str.toString();
        this.splittedSnps = splittedSnps;
    }

    @Override
    public String toString() {
        return minors + snps;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Clique clique = (Clique) o;

        if (minors != null ? !minors.equals(clique.minors) : clique.minors != null) return false;
        if (snps != null ? !snps.equals(clique.snps) : clique.snps != null) return false;
        return splittedSnps != null ? splittedSnps.equals(clique.splittedSnps) : clique.splittedSnps == null;
    }

    @Override
    public int hashCode() {
        int result = minors != null ? minors.hashCode() : 0;
        result = 31 * result + (snps != null ? snps.hashCode() : 0);
        result = 31 * result + (splittedSnps != null ? splittedSnps.hashCode() : 0);
        return result;
    }
}
