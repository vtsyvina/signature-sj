package by.bsu.algorithms;

import by.bsu.distance.HammingDistance;
import by.bsu.model.Clique;
import by.bsu.model.IlluminaSNVSample;
import by.bsu.model.PairEndRead;
import by.bsu.model.SNVResultContainer;
import by.bsu.model.SNVStructure;
import by.bsu.start.Start;
import by.bsu.util.AlgorithmUtils;
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

public class SNVIlluminaMethod extends AbstractSNV {
    private static String al = "ACGT-N";
    public static int minorCount = al.length() - 2;
    public static List<Clique> savageHaplo;

    /**
     * Main method for SNV method on PacBio reads input
     *
     * @param sample Given input reads with N in the beginning or end of read to make them all equal size see in {@link DataReader#readOneLined}
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public Set<SNVResultContainer> getHaplotypes(IlluminaSNVSample sample) {
        long start = System.currentTimeMillis();
        System.out.println("Start 2SNV method");
        System.out.print("Compute profile");
        double[][] profile = Utils.profile(sample, al);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute split columns");
        IlluminaSNVSample splitted = DataReader.splitColumns(sample, profile, "ACGT-");
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute SNV data structure");
        SNVStructure structure = SNVStructureBuilder.buildIllumina(splitted, sample, profile);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute cliques");
        String consensus = Utils.consensus(structure.profile, al);
        //TODO remove it
        savageHaplo = new ArrayList<>();
        for (String s : Start.savage) {
            Set<Integer> snps = new HashSet<>();
            for (int i = 0; i < Math.min(s.length(), consensus.length()); i++) {
                if (s.charAt(i) != consensus.charAt(i)) {
                    snps.add(splittedPosition(i, s.charAt(i), structure));
                }
            }
            savageHaplo.add(new Clique(snps, structure));
        }
        //getMergedCluques first time to get cliques
        Set<Set<Integer>> cliques = run(splitted, structure, sample, true);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Remove reads with low quality ");
        //remove bad reads( >23 mistakes outside of cliques positions)
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / 4)).distinct().sorted().collect(Collectors.toList());
        List<PairEndRead>[] mistakes = new ArrayList[splitted.referenceLength];
        HammingDistance h = new HammingDistance();
        for (PairEndRead read : sample.reads) {
            int apply = h.apply(read.l, consensus.substring(read.lOffset, read.lOffset + read.l.length()));
            if (read.r.length() > 0) {
                apply += h.apply(read.r, consensus.substring(read.rOffset, read.rOffset + read.r.length()));
            }
            for (int i = 0; i < read.l.length(); i++) {
                if (read.l.charAt(i) == 'N') {
                    apply--;
                }
            }
            for (int i = 0; i < read.r.length(); i++) {
                if (read.r.charAt(i) == 'N') {
                    apply--;
                }
            }
            for (Integer i : allPositionsInCliques) {
                if (readCharAtPosition(read, i) != 'N' && readCharAtPosition(read, i) != consensus.charAt(i)) {
                    apply--;
                }
            }
            if (mistakes[apply] == null) {
                mistakes[apply] = new ArrayList<>();
            }
            mistakes[apply].add(read);
        }
        int readsNumber = sample.reads.size();
        sample.reads = new ArrayList<>();
        for (int i = 0; i < 5; i++) {
            if (mistakes[i] != null) {
                sample.reads.addAll(mistakes[i]);
            }
        }
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute profile");
        //after removing all bad reads, rerun whole process again
        profile = Utils.profile(sample, al);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute split columns");
        splitted = DataReader.splitColumns(sample, profile, "ACGT-");
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute SNV data structure");
        structure = SNVStructureBuilder.buildIllumina(splitted, sample, profile);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute cliques");
        cliques = run(splitted, structure, sample, false);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Start getting haplotypes");
        // divide by clusters and find haplotypes
        Set<SNVResultContainer> snvResultContainers = processCliques(cliques, structure, sample, false);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));


        return snvResultContainers;
    }


    public Set<Set<Integer>> run(IlluminaSNVSample splittedSample, SNVStructure struct, IlluminaSNVSample src, boolean log) {
        int count = 0;
        int more = 0;
        String consensus = Utils.consensus(struct.profile, al);
        List<Set<Integer>> adjacencyList = new ArrayList<>();
        //maybe we will remove second minors
        int[] firstMinors = new int[src.referenceLength];
        for (int i = 0; i < firstMinors.length; i++) {
            double max = 0;
            double minor = 0;
            int maxJ = 0;
            int minorJ = 0;
            for (int j = 0; j < struct.profile.length; j++) {
                if (struct.profile[j][i] > max) {
                    minorJ = maxJ;
                    maxJ = j;
                    minor = max;
                    max = struct.profile[j][i];
                } else if (struct.profile[j][i] > minor) {
                    minorJ = j;
                    minor = struct.profile[j][i];
                }
            }
            firstMinors[i] = minorJ;
        }
        for (int i = 0; i < splittedSample.referenceLength; i++) {
            adjacencyList.add(new HashSet<>());
            int l = struct.rowMinors[i].length;
            if (l < 10) {
                continue;
            }
            int[] hits = getHits(struct, struct.rowMinors[i], splittedSample.referenceLength);
            for (int j = 0; j < hits.length; j++) {
                //skip small amount of hits
                if (i != j && hits[j] >= 10 && Math.abs(i - j) >= 20) {
                    //get unsplitted columns, minors, o_kl
                    int first = i / minorCount;
                    int second = j / minorCount;
                    int allele1 = getAllele(i, struct);
                    int allele2 = getAllele(j, struct);

                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    //skip second and other less frequent minors
//                    if (allele1 != firstMinors[first] || allele2 != firstMinors[second]){
//                        continue;
//                    }
                    if (m1 == '-' || m2 == '-') {
                        continue;
                    }
                    double hitsM = 100 * hits[j] / Math.min((double) l, (double) struct.rowMinors[j].length);
                    /*
                     * false 1 means that in actual sample it has another minor or N in given position
                     */
                    int o22 = hits[j];
                    int o21 = struct.rowMinors[i].length; //all 2*
                    int o12 = struct.rowMinors[j].length; //all *2
                    // subtract 2N and false 21 from o21
                    for (int k = 0; k < struct.rowMinors[i].length; k++) {
                        PairEndRead read = src.reads.get(struct.rowMinors[i][k]);
                        //find in what part of read is second
                        if (read.lOffset <= second && read.lOffset + read.l.length() > second) {
                            if (read.l.charAt(second - read.lOffset) != consensus.charAt(second)) {
                                o21--;
                            }
                        } else if (read.rOffset <= second && read.rOffset + read.r.length() > second) {
                            if (read.r.charAt(second - read.rOffset) != consensus.charAt(second)) {
                                o21--;
                            }
                        } else {
                            o21--;//second doesn't have this position at all
                        }

                    }
                    if (o21 == 0) {
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e zero",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, 0.0));
                        adjacencyList.get(i).add(j);
                        count++;
                        continue;
                    }
                    //subtract N2 and false 12 from o12
                    for (int k = 0; k < struct.rowMinors[j].length; k++) {
                        PairEndRead read = src.reads.get(struct.rowMinors[j][k]);
                        //find in what part of read is first
                        if (read.lOffset <= first && read.lOffset + read.l.length() > first) {
                            if (read.l.charAt(first - read.lOffset) != consensus.charAt(first)) {
                                o12--;
                            }
                        } else if (read.rOffset <= first && read.rOffset + read.r.length() > first) {
                            if (read.r.charAt(first - read.rOffset) != consensus.charAt(first)) {
                                o12--;
                            }
                        } else {
                            o12--;//first doesn't have this position at all
                        }
                    }
                    if (o12 == 0) {
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e zero",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, 0.0));
                        adjacencyList.get(i).add(j);
                        count++;
                        continue;
                    }

                    //amount of common reads for i and j column (first two summands to get amount of 1 in j(major and all minors that turned to 1)), then add 21 and 22. We need inly to subtract 1N
                    //subtract 1N from reads
                    int reads = getCommonReadsCount(struct.readsAtPosition[first], struct.readsAtPosition[second]);

                    int o11 = getO11(struct, src, consensus, i, j, first, second, o22, o21, o12, reads);
                    if (o11 == 0) {
                        int tmp = o11;
                        o11 = o22;
                        o22 = tmp;
                    }
                    //start calculate p-value, starting with p
                    //double p = struct.rowMinors[i].length/(double)(sample.reads.length - struct.rowN[first].length);
                    double p = (o12 * o21) / ((double) o11 * reads);

                    if (p > 1) {
                        more++;
                        continue;
                    }
                    if (p < 1E-12) {
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e zero",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, p));
                        adjacencyList.get(i).add(j);
                        count++;
                    } else {
                        double pvalue = Utils.binomialPvalue(o22, p, reads);
                        if (pvalue < 0.0000001 / (splittedSample.referenceLength * (splittedSample.referenceLength - 1) / 2)) {
                            if (log)
                                System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%.2f p=%.3e",
                                        first, second, m1, m2, l, struct.rowMinors[j].length, hitsM, pvalue));
                            adjacencyList.get(i).add(j);
                            count++;
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
        if (log) System.out.println(adjacencyList.stream().mapToInt(Set::size).sum());
        if (log) System.out.println(more);
        Set<Set<Integer>> mergedCliques = getMergedCliques(adjacencyList, struct);
        //Set<Set<Integer>> mergedCliques1 = super.getMergedCliques(adjacencyList);
        int conflictsResolved = 0;
        //resolve conflicts in cliques, leave only snps on the same position with higher frequency
        for (Set<Integer> clique : mergedCliques) {
            Iterator<Integer> it = clique.iterator();
            while (it.hasNext()) {
                Integer v = it.next();
                int row = v - v % 4;
                int allele = getAllele(v, struct);
                for (int i = 0; i < minorCount; i++) {
                    if (row + i == v) {
                        continue;
                    }
                    if (clique.contains(row + i)) {
                        if (struct.profile[allele][v / 4] < struct.profile[getAllele(row + i, struct)][v / 4]) {
                            it.remove();
                            conflictsResolved++;
                        }
                    }
                }
            }
        }
        System.out.println("Conflicts in cliques resolved " + conflictsResolved);
        //List<Clique> collect = mergedCliques.stream().map(c -> new Clique(c, struct)).sorted(Comparator.comparingInt(c -> c.minors.length())).collect(Collectors.toList());
        return mergedCliques;
    }


    protected Set<Set<Integer>> getMergedCliques(List<Set<Integer>> adjacencyList, SNVStructure structure) {
        int edges = 0;
        Set<Integer> p = new HashSet<>();
        for (int i = 0; i < adjacencyList.size(); i++) {
            p.add(i);
        }
        List<Set<Integer>> cliques = AlgorithmUtils.findCliques(new HashSet<>(), p, new HashSet<>(), adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toList());
        System.out.println("Cliques before merge " + cliques.size());
        List<Set<Integer>> cliquesAdjacencyMatrix = new ArrayList<>();
        cliques.forEach(i -> cliquesAdjacencyMatrix.add(new HashSet<>()));
        for (int i = 0; i < cliques.size(); i++) {
            for (int j = 0; j < cliques.size(); j++) {
                if (i == j) {
                    continue;
                }
                if (mergeScore(cliques.get(i), cliques.get(j), adjacencyList) > 2) {
                    cliquesAdjacencyMatrix.get(i).add(j);
                    cliquesAdjacencyMatrix.get(j).add(i);
                    edges++;
                }
            }
        }
        int[] components = AlgorithmUtils.connectedComponents(cliquesAdjacencyMatrix);
        Map<Integer, Set<Integer>> mergedCliques = new HashMap<>();
        for (int i = 0; i < components.length; i++) {
            if (!mergedCliques.containsKey(components[i])) {
                mergedCliques.put(components[i], cliques.get(i));
            } else {
                mergedCliques.get(components[i]).addAll(cliques.get(i));
            }
        }
        System.out.println("Edges between cliques in total " + edges);
        System.out.println("Merged cliques " + mergedCliques.size());
        return new HashSet<>(mergedCliques.values());
    }

    private int getO11(SNVStructure struct, IlluminaSNVSample src, String consensus, int i, int j, int first, int second, int o22, int o21, int o12, int reads) {
        int o11 = reads - o12 - o21 - o22; //all 11(and some false positive 11),12
        //1 1'
        for (int k = 0; k < minorCount; k++) {
            if (k == j % minorCount) {
                continue;
            }
            int currentMinor = j - j % 4 + k;
            for (int m = AlgorithmUtils.binarySearch(struct.rowMinors[currentMinor], struct.readsAtPosition[second][0]); m < struct.rowMinors[currentMinor].length; m++) {
                PairEndRead read = src.reads.get(struct.rowMinors[currentMinor][m]);
                //find in what part of read is first
                if (read.lOffset <= first && read.lOffset + read.l.length() > first) {
                    if (read.l.charAt(first - read.lOffset) == consensus.charAt(first)) {
                        o11--;
                    }
                } else if (read.rOffset <= first && read.rOffset + read.r.length() > first) {
                    if (read.r.charAt(first - read.rOffset) == consensus.charAt(first)) {
                        o11--;
                    }
                }
            }

        }
        //1' 1
        for (int k = 0; k < minorCount; k++) {
            if (k == i % minorCount) {
                continue;
            }
            int currentMinor = i - i % 4 + k;
            for (int m = AlgorithmUtils.binarySearch(struct.rowMinors[currentMinor], struct.readsAtPosition[first][0]); m < struct.rowMinors[currentMinor].length; m++) {
                PairEndRead read = src.reads.get(struct.rowMinors[currentMinor][m]);
                //find in what part of read is second
                if (read.lOffset <= second && read.lOffset + read.l.length() > second) {
                    if (read.l.charAt(second - read.lOffset) == consensus.charAt(second)) {
                        o11--;
                    }
                } else if (read.rOffset <= second && read.rOffset + read.r.length() > second) {
                    if (read.r.charAt(second - read.rOffset) == consensus.charAt(second)) {
                        o11--;
                    }
                }
            }
        }
        //1' 1', 1' 2, 2 1'
        for (int k = 0; k < minorCount; k++) {
            if (struct.rowMinors[i - i % 4 + k].length == 0) {
                continue;
            }
            for (int l = 0; l < minorCount; l++) {
                if ((l == j % minorCount && k == i % minorCount) || struct.rowMinors[j - j % 4 + l].length == 0) {
                    continue;
                }
                o11 -= getSortedArraysIntersectionCount(struct.rowMinors[i - i % 4 + k], struct.rowMinors[j - j % 4 + l]);
            }
        }
        return o11;
    }

    private int getCommonReadsCount(int[] readsAtFirstPosition, int[] readsAtSecondPosition) {
        return getSortedArraysIntersectionCount(readsAtFirstPosition, readsAtSecondPosition);
    }

    private int getSortedArraysIntersectionCount(int[] x, int[] y) {
        int result = 0;
        int si = AlgorithmUtils.binarySearch(y, x[0]);
        int fi = 0;
        int firstEnd = x.length;
        int secondEnd = y.length;
        while (fi < firstEnd && si < secondEnd) {
            if (x[fi] < y[si]) {
                fi++;
            } else if (x[fi] > y[si]) {
                si++;
            } else {
                result++;
                fi++;
                si++;
            }
        }
        return result;
    }


    /**
     * Check that given read has major allele on given position
     *
     * @param consensus Whole sample consensus
     * @param position  Position to check
     * @param o11       current o11 value to correct
     * @param read      given read
     * @return o11 if read doesn't have major on given position, o11-1 otherwise
     */
    private int checkTrueMajor(String consensus, int position, int o11, PairEndRead read) {
        if (read.lOffset <= position && read.lOffset + read.l.length() > position && read.l.charAt(position - read.lOffset) == consensus.charAt(position)) {
            o11--; //subtract only if read has i column and major on it
        } else if (read.rOffset <= position && read.rOffset + read.r.length() > position && read.r.charAt(position - read.rOffset) == consensus.charAt(position)) {
            o11--;
        }
        return o11;
    }

    /**
     * If left and right reads overlap - merge them into just left read
     * @param reads
     * @return
     */
    public List<PairEndRead> processOverlaps(List<PairEndRead> reads) {
        for (PairEndRead read : reads) {
            if (read.lOffset + read.l.length() > read.rOffset && read.rOffset != -1) {
                StringBuilder newL = new StringBuilder(read.l);
                //maybe r lies entirely in l, then no copying
                if (read.lOffset + read.l.length() - read.rOffset < read.r.length()) {
                    newL.append(read.r.substring(read.lOffset + read.l.length() - read.rOffset));
                }
                read.l = newL.toString();
                read.r = "";
                read.rOffset = -1;
            }
        }
        return reads;
    }

    private boolean readHasPosition(PairEndRead read, int position) {
        if (read.lOffset <= position && read.lOffset + read.l.length() > position) {
            return true;
        } else if (read.rOffset <= position && read.rOffset + read.r.length() > position) {
            return true;
        }
        return false;
    }

    private boolean readHasCharAtPosition(PairEndRead read, int position, char m) {
        if (read.lOffset <= position && read.lOffset + read.l.length() > position) {
            return read.l.charAt(position - read.lOffset) == m;
        } else if (read.rOffset <= position && read.rOffset + read.r.length() > position) {
            return read.r.charAt(position - read.rOffset) == m;
        }
        return false;
    }

    private char readCharAtPosition(PairEndRead read, int position) {
        if (read.lOffset <= position && read.lOffset + read.l.length() > position) {
            return read.l.charAt(position - read.lOffset);
        } else if (read.rOffset <= position && read.rOffset + read.r.length() > position) {
            return read.r.charAt(position - read.rOffset);
        }
        return 'N';
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
    public Set<SNVResultContainer> processCliques(Set<Set<Integer>> cliques, SNVStructure struct, IlluminaSNVSample src, boolean log) {
        String consensus = Utils.consensus(struct.profile, al);
        cliques.add(new HashSet<>());
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / minorCount)).distinct().sorted().collect(Collectors.toList());
        List<String> allCliquesCharacters = getAllCliquesCharacters(struct, consensus, cliques, allPositionsInCliques);

        Set<Clique> cliquesSet = new HashSet<>();
        cliques.forEach(c -> cliquesSet.add(new Clique(c, struct)));
        Map<String, Set<PairEndRead>> clusters = buildClusters(src, allPositionsInCliques, allCliquesCharacters, consensus);
        //skip clusters with less than 10 reads. Do some stuff for transforming output into human-friendly format
        Set<SNVResultContainer> haplotypes = clusters.entrySet().stream().filter(s -> s.getValue().size() > 10).map(s -> {
            Set<PairEndRead> cluster = s.getValue();
            IlluminaSNVSample snvSample = new IlluminaSNVSample("tmp", new ArrayList<>(cluster), src.referenceLength);
            //if there is no reads, put consensus there
            double[][] profile = Utils.profile(snvSample, al);
            for (int i = 0; i < profile[0].length; i++) {
                double max = 0;
                for (int j = 0; j < profile.length; j++) {
                    if (profile[j][i] > max) {
                        max = profile[j][i];
                    }
                }
                if (!(max > 0)) {
                    int i1 = al.indexOf(consensus.charAt(i)) == -1 ? '-' : al.indexOf(consensus.charAt(i));
                    profile[i1][i] = 1;
                }
            }
            String haplotype = Utils.consensus(profile, al);
            Set<Integer> snps = new HashSet<>();
            for (int i = 0; i < haplotype.length(); i++) {
                if (haplotype.charAt(i) != consensus.charAt(i)) {
                    snps.add(splittedPosition(i, haplotype.charAt(i), struct));
                }
            }
            Clique haplotypeClique = new Clique(snps, struct);
            SNVResultContainer container = new SNVResultContainer(s.getKey(), haplotypeClique, haplotype, cluster);
            container.sourceClique = getSourceClique(allPositionsInCliques, cliquesSet, s.getKey(), container);
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
    private Map<String, Set<PairEndRead>> buildClusters(IlluminaSNVSample src, List<Integer> allPositionsInCliques, List<String> allCliquesCharacters, String consensus) {
        Map<String, Set<PairEndRead>> clusters = new HashMap<>();
        allCliquesCharacters.forEach(s -> clusters.put(s, new HashSet<>()));
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
        for (PairEndRead read : src.reads) {
            List<Integer> distancesFromCliques = new ArrayList<>();
            for (String c : allCliquesCharacters) {
                int d = 0;
                int coincidences = 0;
                boolean isConsensusClique = c.equals(consensusClique);
                for (int i = 0; i < allPositionsInCliques.size(); i++) {
                    char charAtPosition = readCharAtPosition(read, allPositionsInCliques.get(i));
                    //increase distance
                    if (charAtPosition != 'N'){
                        if (charAtPosition != c.charAt(i)) {
                            d++;
                        }
                        //coincidence with current clique
                        if (c.charAt(i) != consensus.charAt(allPositionsInCliques.get(i))) {
                            coincidences++;
                        }
                        //coincidence with any cluque snp. Only for consensus clique
                        if (isConsensusClique){
                            coincidences++;
                        }
                    }
                }

                if (coincidences == 0) {
                    d = 1_000_000;
                }
                distancesFromCliques.add(d);
            }
            int minDistance = distancesFromCliques.stream().mapToInt(c -> c).min().getAsInt();
            if (minDistance < 1_000_000) {
                for (int i = 0; i < distancesFromCliques.size(); i++) {
                    if (distancesFromCliques.get(i) == minDistance) {
                        clusters.get(allCliquesCharacters.get(i)).add(read);
                    }
                }
            }

        }
        return clusters;
    }
}
