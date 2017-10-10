package by.bsu.algorithms;

import by.bsu.distance.HammingDistance;
import by.bsu.model.Clique;
import by.bsu.model.SNVResultContainer;
import by.bsu.model.SNVStructure;
import by.bsu.model.Sample;
import by.bsu.util.DataReader;
import by.bsu.util.Utils;
import by.bsu.util.builders.SNVStructureBuilder;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static by.bsu.util.Utils.distinctByKey;

public class SNVPacBioMethod {
    public static String al = "ACGT-N";
    public static int minorCount = al.length() - 2;

    /**
     * Main method for SNV method on PacBio reads input
     *
     * @param sample Given input reads with N in the beginning or end of read to make them all equal size see in {@link DataReader#readOneLined}
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public Set<SNVResultContainer> getHaplotypes(Sample sample) {
        String consensus = Utils.consensus(sample.sequences, al);
        double[][] profile = Utils.profile(sample, al);
        Sample splitted = DataReader.splitColumns(sample, profile, "ACGT-");
        SNVStructure structure = SNVStructureBuilder.buildPacBio(splitted, sample, profile);
        //run first time to get cliques
        Set<Set<Integer>> cliques = run(splitted, structure, sample, false);

        //remove bad reads( >23 mistakes outside of cliques positions)
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / 4)).distinct().sorted().collect(Collectors.toList());
        int[] mistakes = new int[2000];
        List<String> newSequences = new ArrayList<>();
        for (String sequence : sample.sequences) {
            int apply = new HammingDistance().apply(sequence, consensus);
            for (int i = 0; i < sequence.length(); i++) {
                if (sequence.charAt(i) == 'N') {
                    apply--;
                }
            }
            for (Integer i : allPositionsInCliques) {
                if (sequence.charAt(i) != 'N' && sequence.charAt(i) != consensus.charAt(i)) {
                    apply--;
                }
            }
            mistakes[apply]++;
            //TODO remove only 10 percent
            if (apply <= 23) {
                newSequences.add(sequence);
            }
        }
        sample.sequences = newSequences.toArray(new String[newSequences.size()]);
        //after removing all bad reads, rerun whole process again
        profile = Utils.profile(sample, al);
        splitted = DataReader.splitColumns(sample, profile, "ACGT-");
        structure = SNVStructureBuilder.buildPacBio(splitted, sample, profile);
        cliques = run(splitted, structure, sample, true);
        // divide by clusters and find haplotypes
        return processCliques(cliques, structure, sample, true);
    }

    /**
     * Auxiliary method that calculates all hits(and p-value) between all alleles and create cliques according to them
     * It will merge cliques that share more than 50% of edges
     *
     * @param splittedSample splitted sample where each row substitutes with 'minorCount'(4 in our case) where corresponding minor is denoted as 2 and all other alleles as 1
     * @param struct         SNV structure that helps to count hits between minors {@link SNVStructureBuilder} {@link SNVStructure}
     * @param src            Source sample with given PacBio reads
     * @param log            boolean value if we want to see some additional info during algorithm's work (debug purposes)
     * @return Set of SNPs position for all found cliques
     */
    public Set<Set<Integer>> run(Sample splittedSample, SNVStructure struct, Sample src, boolean log) {
        //TODO remove allHits
        int[][] allHits = new int[splittedSample.sequences[0].length()][];
        String consensus = Utils.consensus(struct.profile, al);
        List<Set<Integer>> adjacencyMatrix = new ArrayList<>();
        for (int i = 0; i < splittedSample.sequences[0].length(); i++) {
            adjacencyMatrix.add(new HashSet<>());
            int l = struct.rowMinors[i].length;
            if (l < 10) {
                continue;
            }
            int[] hits = getHits(struct, struct.rowMinors[i], splittedSample.sequences[i].length());
            allHits[i] = hits;
            for (int j = 0; j < hits.length; j++) {
                //skip small amount of hits
                if (i != j && hits[j] >= 10 && Math.abs(i - j) > 20) {
                    //get unsplitted columns, minors, o_kl
                    int first = i / minorCount;
                    int second = j / minorCount;
                    int allele1 = i % minorCount >= Utils.getMajorAllele(struct.profile, first) ?
                            i % minorCount + 1 : i % minorCount;
                    int allele2 = j % minorCount >= Utils.getMajorAllele(struct.profile, second) ?
                            j % minorCount + 1 : j % minorCount;

                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    if (m1 == '-' || m2 == '-') {
                        continue;
                    }
                    /*
                     * false 1 means that in actual sample it has another minor or N in given position
                     */
                    int o22 = hits[j];
                    int o21 = struct.rowMinors[i].length; //all 2*
                    int o12 = struct.rowMinors[j].length; //all *2
                    // subtract all that a not 21 from o21
                    for (int k = 0; k < struct.rowMinors[i].length; k++) {
                        if (src.sequences[struct.rowMinors[i][k]].charAt(second) != consensus.charAt(second)) {
                            o21--;
                        }
                    }
                    // subtract all that a not 12 from o12
                    for (int k = 0; k < struct.rowMinors[j].length; k++) {
                        if (src.sequences[struct.rowMinors[j][k]].charAt(first) != consensus.charAt(first)) {
                            o12--;
                        }
                    }

                    int o11 = struct.majorsInRow[first]; //all 11(and some false positive 11),12,1N
                    //subtract 1N from o11
                    for (int k = 0; k < struct.rowN[second].length; k++) {
                        if (src.sequences[struct.rowN[second][k]].charAt(first) == consensus.charAt(second)) {//make sure that 1 in minColumn is true 1
                            o11--;
                        }
                    }
                    //subtract 12 and false positive 11 from o11
                    for (int k = 0; k < minorCount; k++) {
                        int[] rowMinor = struct.rowMinors[second * minorCount + k];
                        for (int m : rowMinor) {
                            if (src.sequences[m].charAt(first) == consensus.charAt(first)) { //make sure that 1 in minColumn is true 1
                                o11--;
                            }
                        }
                    }
                    //amount of common reads for i and j column
                    int reads = src.sequences[0].length();
                    //start calculate p-value, starting with p
                    //double p = struct.rowMinors[i].length/(double)(sample.reads.length - struct.rowN[first].length);
                    double p = (o12 * o21) / ((double) o11 * reads);
                    //p = struct.profile[allele1][first] * struct.profile[allele2][second];
                    if (p < 1E-12) {
                        double hitsM = 100 * hits[j] / Math.min((double) l, (double) struct.rowMinors[j].length);
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e zero",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, p));
                        adjacencyMatrix.get(i).add(j);
                    } else {
                        double pvalue = Utils.binomialPvalue(o22, p, reads);
                        double hitsM = 100 * hits[j] / Math.min((double) l, (double) struct.rowMinors[j].length);
                        if (pvalue < 0.00001 / (reads * (reads - 1) / 2)) {
                            if (log)
                                System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e",
                                        first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, pvalue));
                            adjacencyMatrix.get(i).add(j);
                        }
                    }

                }
            }
        }

        /*
         *   if for any 2 positions we have edges like  X <-> Y and X <-> Z,
         *   then we delete edge with less frequency of second allele(to avoid false positive cliques)
         */
        for (int i = 0; i < splittedSample.sequences[0].length(); i++) {
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
        return getMergedCliques(adjacencyMatrix);

    }

    /**
     * Methods process computed cliques. At first it separate reads into clusters and then compute consensus for each cluster and some additional information
     *
     * @param cliques Set of cliques in form of splitted SNPs(where position and minor are encoded in a single number)
     * @param struct  SNV structure
     * @param src     Source sample
     * @param log     boolean value if we want to see some additional info during algorithm's work (debug purposes)
     * @return Set of containers with results. Each container contains haplotype itself and some additional helpful information
     */
    public Set<SNVResultContainer> processCliques(Set<Set<Integer>> cliques, SNVStructure struct, Sample src, boolean log) {
        String consensus = Utils.consensus(struct.profile, al);
        cliques.add(new HashSet<>());
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / minorCount)).distinct().sorted().collect(Collectors.toList());
        List<String> allCliquesCharacters = getAllCliquesCharacters(struct, consensus, cliques, allPositionsInCliques);

        Set<Clique> cliquesSet = new HashSet<>();
        cliques.forEach(c -> cliquesSet.add(new Clique(c, struct)));
        Map<String, Set<String>> clusters = buildClusters(src, allPositionsInCliques, allCliquesCharacters);
        //skip clusters with less than 10 reads. Do some stuff for transforming output into human-friendly format
        Set<SNVResultContainer> haplotypes = clusters.entrySet().stream().filter(s -> s.getValue().size() > 10).map(s -> {
            Set<String> cluster = s.getValue();
            String haplotype = Utils.consensus(cluster.toArray(new String[cluster.size()]), al);
            Set<Integer> snps = new HashSet<>();
            for (int i = 0; i < haplotype.length(); i++) {
                if (haplotype.charAt(i) != consensus.charAt(i)) {
                    snps.add(splittedPosition(i, haplotype.charAt(i), struct));
                }
            }
            Clique haplotypeClique = new Clique(snps, struct);
            SNVResultContainer container = new SNVResultContainer(s.getKey(), cluster, haplotypeClique, haplotype);
            for (Clique clique : cliquesSet) {
                boolean fl = true;
                for (int i = 0; i < clique.snps.size(); i++) {
                    if (s.getKey().charAt(allPositionsInCliques.indexOf(clique.snps.get(i))) != clique.minors.charAt(i)) {
                        fl = false;
                        break;
                    }
                }
                if (fl) {
                    container.sourceClique = clique;
                    //small hack. If clique is empty than fl will be true for any source clique
                    if (clique.minors.length() > 0) {
                        break;
                    }
                }
            }
            return container;
        }).collect(Collectors.toSet());
        return haplotypes.stream().filter(distinctByKey(p -> p.haplotype)).collect(Collectors.toSet());
    }

    /**
     * Build read clusters based on cliques. Each read will go to nearest clique in terms of Hamming distance
     *
     * @param src                   Source sample
     * @param allPositionsInCliques sorted array or all positions with at least one clique
     * @param allCliquesCharacters  characters in cliques according to allPositionsInCliques.
     *                              Has consensus allele if clique doesn't include particular position from allPositionsInCliques
     * @return Map with clusters, where key is string of clique characters, value is a set of reads
     */
    private Map<String, Set<String>> buildClusters(Sample src, List<Integer> allPositionsInCliques, List<String> allCliquesCharacters) {
        Map<String, Set<String>> clusters = new HashMap<>();
        allCliquesCharacters.forEach(s -> clusters.put(s, new HashSet<>()));
        class DistanceContainer {
            public String sequence;
            public int distance;

            private DistanceContainer(String key, int distance) {
                this.sequence = key;
                this.distance = distance;
            }
        }
        for (String s : src.sequences) {
            List<DistanceContainer> distancesFromCliques = new ArrayList<>();
            for (String c : allCliquesCharacters) {
                int d = 0;
                for (int i = 0; i < allPositionsInCliques.size(); i++) {
                    if (s.charAt(allPositionsInCliques.get(i)) != c.charAt(i)) {
                        d++;
                    }
                }
                distancesFromCliques.add(new DistanceContainer(s, d));
            }
            int minDistance = distancesFromCliques.stream().mapToInt(c -> c.distance).min().getAsInt();
            for (int i = 0; i < distancesFromCliques.size(); i++) {
                DistanceContainer c = distancesFromCliques.get(i);
                if (c.distance == minDistance) {
                    clusters.get(allCliquesCharacters.get(i)).add(c.sequence);
                }
            }
        }
        return clusters;
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
    private List<String> getAllCliquesCharacters(SNVStructure struct, String consensus, Set<Set<Integer>> cliques, List<Integer> allPositionsInCliques) {
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

    /**
     * Method to build all cliques based on adjacencyMatrix for SNPs.
     * Also it merges cliques if they have more than 50% of edges in common into 'pseudoclique'
     *
     * @param adjacencyMatrix Matrix with edges (in splitted sample)
     * @return set of cliques (clique - set of splitted positions that encode position + minor)
     */
    private Set<Set<Integer>> getMergedCliques(List<Set<Integer>> adjacencyMatrix) {
        Set<Integer> p = new HashSet<>();
        for (int i = 0; i < adjacencyMatrix.size(); i++) {
            p.add(i);
        }
        Set<Set<Integer>> cliques = new BronKerbosch().findCliques(new HashSet<>(), p, new HashSet<>(), adjacencyMatrix).stream().filter(c -> c.size() > 1).collect(Collectors.toSet());
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
                    int score = mergeScore(clique, clique2, adjacencyMatrix);
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
     * Calculates hits for given row of minors.
     * Gets all minors in column for a certain position for minorRow and increments appropriate hits position
     *
     * @param struct          SNV data structure
     * @param rowMinor        array for all minors for some position
     * @param referenceLength Reference length
     * @return array with all minor hits for given row
     */
    private int[] getHits(SNVStructure struct, int[] rowMinor, int referenceLength) {
        int[] hits = new int[referenceLength];

        for (int j = 0; j < rowMinor.length; j++) {
            int[] column = struct.colMinors[rowMinor[j]];
            for (int aColumn : column) {
                hits[aColumn]++;
            }
        }
        return hits;
    }

    //TODO remove method, it's just for debugging
    private int getHitsBetweenPositions(int i, int j, char m1, char m2, int[][] allHits, SNVStructure struct, Sample src) {
        int allele1 = al.indexOf(m1) >= Utils.getMajorAllele(struct.profile, i) ? al.indexOf(m1) - 1 : al.indexOf(m1);
        int allele2 = al.indexOf(m2) >= Utils.getMajorAllele(struct.profile, j) ? al.indexOf(m2) - 1 : al.indexOf(m2);
        int first = i * minorCount + allele1;
        int second = j * minorCount + allele2;
        System.out.println("positions of hits between " + i + " " + m1 + " " + j + " " + m2 + ":");
        for (int k = 0; k < struct.rowMinors[first].length; k++) {
            if (src.sequences[struct.rowMinors[first][k]].charAt(j) == m2) {
                System.out.println(struct.rowMinors[first][k]);
            }
        }

        return allHits[first][second];
    }

    /**
     * Gives corresponding minot for given splitted snp
     */
    private char minor(int splittedSnp, SNVStructure structure) {
        int allele1 = splittedSnp % minorCount >= Utils.getMajorAllele(structure.profile, splittedSnp / minorCount) ?
                splittedSnp % minorCount + 1 :
                splittedSnp % minorCount;
        return al.charAt(allele1);
    }

    /**
     * Transforms position + minor into splitted SNP position
     */
    private int splittedPosition(int pos, char minor, SNVStructure structure) {
        int r = pos * minorCount;
        int allele = al.indexOf(minor) >= Utils.getMajorAllele(structure.profile, pos) ? al.indexOf(minor) - 1 : al.indexOf(minor);
        return r + allele;
    }

    /**
     * Compute amount of common edges between the two given cliques. Returns 0 if cliques have common position but different minors
     */
    private int mergeScore(Set<Integer> c1, Set<Integer> c2, List<Set<Integer>> adjacencyMatrix) {
        int matches = 0;
        for (Integer i : c1) {
            for (Integer i2 : c2) {
                if (i / minorCount == i2 / minorCount && i % minorCount != i2 % minorCount) {
                    return 0;
                }
                if (adjacencyMatrix.get(i).contains(i2)) {
                    matches++;
                }
            }
        }
        return matches;
    }

    /**
     * Simple implementation of Bron - Kerbosch algorithm for finding all maximum cliques (psedocode is on Wiki)
     */
    class BronKerbosch {
        Set<Set<Integer>> findCliques(Set<Integer> clique, Set<Integer> p, Set<Integer> x, List<Set<Integer>> adjacencyMatrix) {
            Set<Set<Integer>> result = new HashSet<>();
            if (p.isEmpty() && x.isEmpty()) {
                result.add(clique);
            }
            Iterator<Integer> it = p.iterator();
            while (it.hasNext()) {
                Integer v = it.next();
                Set<Integer> newCLique = new HashSet<>(clique);
                newCLique.add(v);
                Set<Integer> pIntersection = new HashSet<>(p);
                Set<Integer> xIntersection = new HashSet<>(x);
                pIntersection.retainAll(adjacencyMatrix.get(v));
                xIntersection.retainAll(adjacencyMatrix.get(v));
                result.addAll(findCliques(newCLique, pIntersection, xIntersection, adjacencyMatrix));
                it.remove();
                x.add(v);
            }
            return result;
        }
    }
}


