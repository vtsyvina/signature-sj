package by.bsu.algorithms;

import by.bsu.algorithms.EM.PacBioEM;
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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static by.bsu.util.Utils.distinctByKey;

public class SNVPacBioMethod extends AbstractSNV {


    /**
     * Main method for SNV method on PacBio reads input
     *
     * @param sample Given input reads with N in the beginning or end of read to make them all equal size see in {@link DataReader#readOneLined}
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public List<SNVResultContainer> getHaplotypes(Sample sample, boolean log) {
        String consensus = Utils.consensus(sample.sequences, al);
        if (log) System.out.println("Start 2SNV method");
        if (log) System.out.print("Compute profile");
        double[][] profile = Utils.profile(sample, al);
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute split columns");
        Sample splitted = DataReader.splitColumns(sample, profile, "ACGT-");
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute SNV data structure");
        SNVStructure structure = SNVStructureBuilder.buildPacBio(splitted, sample, profile);
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute cliques");
        //run first time to get cliques
        Set<Set<Integer>> cliques = run(splitted, structure, sample, log);
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Remove reads with low quality ");
        //remove bad reads( >23 mistakes outside of cliques positions)
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / 4)).distinct().sorted().collect(Collectors.toList());
        int[] mistakes = new int[2000];
        List<String>[] m = new ArrayList[2000];
        for (int i = 0; i < 2000; i++) {
            m[i] = new ArrayList<>();
        }
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
            m[apply].add(sequence);
        }
        outer: for (List<String> strings : m) {
            for (String r : strings) {
                if (newSequences.size() < 0.9*sample.sequences.length){
                    newSequences.add(r);
                }else {
                    break outer;
                }
            }
        }
        Sample newSample = new Sample(sample.name, newSequences.toArray(new String[newSequences.size()]));
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute profile");
        //after removing all bad reads, rerun whole process again
        profile = Utils.profile(newSample, al);
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute split columns");
        splitted = DataReader.splitColumns(newSample, profile, "ACGT-");
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute SNV data structure");
        SNVStructure newStructure = SNVStructureBuilder.buildPacBio(splitted, newSample, profile);
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Compute cliques");
        cliques = run(splitted, newStructure, newSample, log);
        if (log) System.out.println(" - DONE");
        if (log) System.out.print("Start getting haplotypes");
        // divide by clusters and find haplotypes
        List<SNVResultContainer> snvResultContainers = processCliques(cliques, structure, sample, log);
        if (log) System.out.println(" - DONE");
        return snvResultContainers;
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
        List<Set<Integer>> adjacencyList = new ArrayList<>();
        for (int i = 0; i < splittedSample.sequences[0].length(); i++) {
            adjacencyList.add(new HashSet<>());
            int l = struct.rowMinors[i].length;
            if (l < 10) {
                continue;
            }
            int[] hits = getHits(struct, struct.rowMinors[i], splittedSample.sequences[0].length());
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
                    //amount of common reads for i and j column (first two summands to get amount of 1 in j(major and all minors that turned to 1)), then add 21 and 22. We need inly to subtract 1N
                    int reads = src.sequences.length - struct.rowN[first].length + o22 + o21;
                    int ni = 0;
                    int nj = 0;
                    int nn = 0;
                    while (ni < struct.rowN[first].length && nj < struct.rowN[second].length) {
                        if (struct.rowN[first][ni] < struct.rowN[second][nj]) {
                            ni++;
                        } else if (struct.rowN[first][ni] < struct.rowN[second][nj]) {
                            nj++;
                        } else {
                            nn++;
                            ni++;
                            nj++;
                        }
                    }
                    //subtract 1N from reads
                    reads -= struct.rowN[second].length - nn;
                    //start calculate p-value, starting with p
                    //double p = struct.rowMinors[i].length/(double)(sample.reads.length - struct.rowN[first].length);
                    double p = (o12 * o21) / ((double) o11 * reads);
                    //p = struct.profile[allele1][first] * struct.profile[allele2][second];
                    double hitsM = 100 * hits[j] / Math.min((double) l, (double) struct.rowMinors[j].length);
                    if (p < 1E-12) {
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e zero",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, p));
                        adjacencyList.get(i).add(j);
                    } else {
                        double pvalue = Utils.binomialPvalue(o22, p, reads);
                        if (pvalue < 0.00001 / (reads * (reads - 1) / 2)) {
                            if (log)
                                System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e",
                                        first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, pvalue));
                            adjacencyList.get(i).add(j);
                        }
                    }

                }
            }
        }

        /*
         *   if for any 2 positions we have edges like  X <-> Y and X <-> Z,
         *   then we delete edge with less frequency of second allele(to avoid false positive cliques)
         */
        removeEdgesForSecondMinors(adjacencyList, struct, log);
        if (log) System.out.println("Edges found "+adjacencyList.stream().mapToInt(Set::size).sum());
        return getMergedCliques(adjacencyList);

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
    public List<SNVResultContainer> processCliques(Set<Set<Integer>> cliques, SNVStructure struct, Sample src, boolean log) {
        String consensus = Utils.consensus(struct.profile, al);
        cliques.add(new HashSet<>());
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / minorCount)).distinct().sorted().collect(Collectors.toList());
        List<String> allCliquesCharacters = getAllCliquesCharacters(struct, consensus, cliques, allPositionsInCliques);

        Set<Clique> cliquesSet = new HashSet<>();
        cliques.forEach(c -> cliquesSet.add(new Clique(c, struct)));
        Map<String, List<String>> clusters = buildClusters(src, allPositionsInCliques, allCliquesCharacters, consensus);
        //skip clusters with less than 10 reads. Do some stuff for transforming output into human-friendly format
        List<SNVResultContainer> haplotypes = clusters.entrySet().stream().filter(s -> s.getValue().size() > 10).map(s -> {
            List<String> cluster = s.getValue();
            String haplotype = Utils.consensus(cluster.toArray(new String[cluster.size()]), al);
            Set<Integer> snps = new HashSet<>();
            for (int i = 0; i < haplotype.length(); i++) {
                if (haplotype.charAt(i) != consensus.charAt(i)) {
                    snps.add(splittedPosition(i, haplotype.charAt(i), struct));
                }
            }
            Clique haplotypeClique = new Clique(snps, struct);
            SNVResultContainer container = new SNVResultContainer(s.getKey(), cluster, haplotypeClique, haplotype);
            container.sourceClique = getSourceClique(allPositionsInCliques, cliquesSet, s.getKey(), container);
            return container;
        }).collect(Collectors.toList());
        List<SNVResultContainer> result = haplotypes.stream().filter(distinctByKey(p -> p.haplotype)).collect(Collectors.toList());
        List<String> h = result.stream().map(s -> s.haplotype).collect(Collectors.toList());
        List<Double> frequencies = new PacBioEM().frequencies(h, src);
        for (int i = 0; i < frequencies.size(); i++) {
            result.get(i).frequency = frequencies.get(i);
        }
        return result.stream().sorted((s1,s2) -> -Double.compare(s1.frequency, s2.frequency)).collect(Collectors.toList());
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
    private Map<String, List<String>> buildClusters(Sample src, List<Integer> allPositionsInCliques, List<String> allCliquesCharacters, String consensus) {
        Map<String, List<String>> clusters = new HashMap<>();
        allCliquesCharacters.forEach(s -> clusters.put(s, new ArrayList<>()));
        class DistanceContainer {
            public String sequence;
            public int distance;
            public int coincidences;

            private DistanceContainer(String key, int distance, int coincidences) {
                this.sequence = key;
                this.distance = distance;
                this.coincidences = coincidences;
            }
        }
        String consensusClique = "";
        for (String characters : allCliquesCharacters) {
            boolean fl = true;
            for (int i = 0; i < characters.length(); i++) {
                if (characters.charAt(i) != consensus.charAt(allPositionsInCliques.get(i))) {
                    fl = false;
                    break;
                }
            }
            if (fl) {
                consensusClique = characters;
            }
        }
        for (String s : src.sequences) {
            List<DistanceContainer> distancesFromCliques = new ArrayList<>();
            for (String c : allCliquesCharacters) {
                int d = 0;
                int coincidences = 0;
                boolean isConsensusClique = c.equals(consensusClique);
                for (int i = 0; i < allPositionsInCliques.size(); i++) {
                    if (s.charAt(allPositionsInCliques.get(i)) != c.charAt(i)) {
                        d++;
                    }
                    //coincidence with current clique
                    if (c.charAt(i) != consensus.charAt(allPositionsInCliques.get(i))) {
                        coincidences++;
                    }
                    //coincidence with current clique
                    if (c.charAt(i) != consensus.charAt(allPositionsInCliques.get(i))) {
                        coincidences++;
                    }
                    //coincidence with any cluque snp. Only for consensus clique
                    if (isConsensusClique) {
                        coincidences++;
                    }
                }
                distancesFromCliques.add(new DistanceContainer(s, d, coincidences));
            }
            int minDistance = distancesFromCliques.stream().mapToInt(c -> c.distance).min().getAsInt();
            if (distancesFromCliques.stream().filter(c -> c.distance == minDistance).count() > 1){
                continue;
            }
            for (int i = 0; i < distancesFromCliques.size(); i++) {
                DistanceContainer c = distancesFromCliques.get(i);
                if (c.distance == minDistance && c.distance != c.coincidences) {
                    clusters.get(allCliquesCharacters.get(i)).add(c.sequence);
                }
            }
        }
        return clusters;
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
}


