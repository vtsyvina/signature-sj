package by.bsu.algorithms;

import by.bsu.model.Clique;
import by.bsu.model.PairEndRead;
import by.bsu.model.SNVResultContainer;
import by.bsu.model.SNVStructure;
import by.bsu.util.AlgorithmUtils;
import by.bsu.util.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class AbstractSNV {
    public static String al = "ACGT-N";
    public static int minorCount = al.length() - 2;

    /**
     * Method to build all cliques based on adjacencyList for SNPs.
     * Also it merges cliques if they have more than 50% of edges in common into 'pseudoclique'
     *
     * @param adjacencyList Matrix with edges (in splitted sample)
     * @return set of cliques (clique - set of splitted positions that encode position + minor)
     */
    protected Set<Set<Integer>> getMergedCliques(List<Set<Integer>> adjacencyList) {
        Set<Integer> p = new HashSet<>();
        for (int i = 0; i < adjacencyList.size(); i++) {
            p.add(i);
        }
        Set<Set<Integer>> cliques = AlgorithmUtils.findCliques(new HashSet<>(), p, new HashSet<>(), adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toSet());
        //trying to merge cliques
        Set<Set<Integer>> mergedCliques = new HashSet<>();
        Set<Set<Integer>> alreadyMerged = new HashSet<>();
        boolean fl = true;
        while (fl) {
            fl = false;
            for (Set<Integer> clique : cliques) {
                int bestScore = 0;
                Set<Integer> bestClique = null;
                if (alreadyMerged.contains(clique)) {
                    continue;
                }
                //fing best merging score for current clique
                for (Set<Integer> clique2 : cliques) {
                    if (clique == clique2 || alreadyMerged.contains(clique2)) {
                        continue;
                    }
                    int score = mergeScore(clique, clique2, adjacencyList);
                    if (score > bestScore) {
                        bestScore = score;
                        bestClique = clique2;
                    }
                }
                // merge two cliques if they have more than 50% of edges in common
                if (bestScore > 0 && 2 * bestScore > clique.size() * bestClique.size()) {
                    HashSet<Integer> newClique = new HashSet<>(clique);
                    newClique.addAll(bestClique);
                    mergedCliques.add(newClique);
                    alreadyMerged.add(bestClique);
                    alreadyMerged.add(clique);
                    fl = true;
                } else {
                    mergedCliques.add(clique);
                    alreadyMerged.add(clique);
                }
            }
            cliques = mergedCliques;
            mergedCliques = new HashSet<>();
            alreadyMerged = new HashSet<>();
        }
        return cliques;
    }

    /**
     * Compute amount of common edges between the two given cliques. Returns 0 if cliques have common position but different minors
     */
    protected int mergeScore(Set<Integer> c1, Set<Integer> c2, List<Set<Integer>> adjacencyList) {
        int matches = 0;
        for (Integer i : c1) {
            for (Integer i2 : c2) {
                if (i / minorCount == i2 / minorCount && i % minorCount != i2 % minorCount) {
                    return 0;
                }
                if (adjacencyList.get(i).contains(i2)) {
                    matches++;
                }
            }
        }
        return matches;
    }

    /**
     * Gives corresponding minot for given splitted snp
     */
    protected char minor(int splittedSnp, SNVStructure structure) {
        int allele1 = splittedSnp % minorCount >= Utils.getMajorAllele(structure.profile, splittedSnp / minorCount) ?
                splittedSnp % minorCount + 1 :
                splittedSnp % minorCount;
        return al.charAt(allele1);
    }

    /**
     * Transforms position + minor into splitted SNP position
     */
    protected int splittedPosition(int pos, char minor, SNVStructure structure) {
        int r = pos * minorCount;
        int allele = al.indexOf(minor) >= Utils.getMajorAllele(structure.profile, pos) ? al.indexOf(minor) - 1 : al.indexOf(minor);
        return r + allele;
    }

    /**
     * Calculates hits for given row of minors.
     * Gets all minors in column for a certain position for minorRow and increments appropriate hits position
     *
     * @param struct          SNV data structure
     * @param rowMinor        array for all minors for some position
     * @param referenceLength Reference length
     * @return array with all minor hits for given row
     */
    protected int[] getHits(SNVStructure struct, int[] rowMinor, int referenceLength) {
        int[] hits = new int[referenceLength];

        for (int j = 0; j < rowMinor.length; j++) {
            int[] column = struct.colMinors[rowMinor[j]];
            for (int aColumn : column) {
                hits[aColumn]++;
            }
        }
        return hits;
    }

    /**
     * Transform each clique into a string of equal length. The string corresponds to clique representation in all positions
     * that are covered by any clique. If clique doesn't have some position in it, than consensus allele will be in this position
     *
     * @param struct                SNV structure
     * @param consensus             Consensus string
     * @param cliques               Set of cliques in form of splitted SNPs(where position and minor are encoded in a single number)
     * @param allPositionsInCliques List of all positions covered by any clique
     * @return list with strings representing given cliques
     */
    protected List<String> getAllCliquesCharacters(SNVStructure struct, String consensus, Set<Set<Integer>> cliques, List<Integer> allPositionsInCliques) {
        List<String> allCliquesCharacters = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        allPositionsInCliques.forEach(i -> str.append(consensus.charAt(i)));
        Set<Clique> clicuesSet = new HashSet<>();
        cliques.forEach(c -> clicuesSet.add(new Clique(c, struct)));
        clicuesSet.forEach(c -> {
            StringBuilder tmp = new StringBuilder(str);
            for (int i = 0; i < c.snps.size(); i++) {
                int snpPosition = allPositionsInCliques.indexOf(c.snps.get(i));
                tmp.replace(snpPosition, snpPosition + 1, String.valueOf(c.minors.charAt(i)));
            }
            allCliquesCharacters.add(tmp.toString());
        });
        return allCliquesCharacters;
    }

    /*
     *   if for any 2 positions we have edges like  X <-> Y and X <-> Z,
     *   then we delete edge with less frequency of second allele(to avoid false positive cliques)
     */
    protected void removeEdgesForSecondMinors(List<Set<Integer>> adjacencyMatrix, SNVStructure struct, boolean log) {
        for (int i = 0; i < adjacencyMatrix.size(); i++) {
            //convert adjacencyMatrix for each position into map (position -> all correlated alleles in this position)
            Map<Integer, Set<Integer>> columnsEdges = new HashMap<>();
            adjacencyMatrix.get(i).forEach(j -> {
                if (!columnsEdges.containsKey(j / minorCount)) {
                    columnsEdges.put(j / minorCount, new HashSet<>());
                }
                columnsEdges.get(j / minorCount).add(j);
            });

            int finalI1 = i;
            // for evry position where we have more than 1 edge
            columnsEdges.entrySet().stream().filter(e -> e.getValue().size() > 1).forEach(e -> {
                final double[] max = {0};
                final int[] maxM = {0};
                e.getValue().forEach(m -> {
                    if (struct.profile[al.indexOf(minor(m, struct))][m / minorCount] > max[0]) {
                        max[0] = struct.profile[al.indexOf(minor(m, struct))][m / minorCount];
                        maxM[0] = m;
                    }
                });
                //remove all non-maximum frequency edges
                e.getValue().forEach(m -> {
                    if (m != maxM[0]) {
                        adjacencyMatrix.get(finalI1).remove(m);
                        adjacencyMatrix.get(m).remove(finalI1);
                        if (log)
                            System.out.println("remove " + finalI1 / minorCount + " " + m / minorCount + " " + minor(finalI1, struct) + " " + minor(m, struct));
                    }
                });
            });
        }
    }

    protected int getAllele(int i, SNVStructure struct){
        return i % minorCount >= Utils.getMajorAllele(struct.profile, i /4) ? i % minorCount + 1 : i % minorCount;
    }

    protected Clique getSourceClique(List<Integer> allPositionsInCliques, Set<Clique> cliquesSet, String cliqueString, SNVResultContainer container) {
        Clique sourceClique = container.sourceClique;
        for (Clique clique : cliquesSet) {
            boolean fl = true;
            for (int i = 0; i < clique.snps.size(); i++) {
                if (cliqueString.charAt(allPositionsInCliques.indexOf(clique.snps.get(i))) != clique.minors.charAt(i)) {
                    fl = false;
                    break;
                }
            }
            if (fl) {
                sourceClique = clique;
                //small hack. If clique is empty than fl will be true for any source clique
                if (clique.minors.length() > 0) {
                    break;
                }
            }
        }
        return sourceClique;
    }
}
