package by.bsu.algorithms;

import by.bsu.algorithms.EM.IlluminaEM;
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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import static by.bsu.util.Utils.distinctByKey;

public class SNVIlluminaMethod extends AbstractSNV {
    private static String al = "ACGT-N";
    public static int minorCount = al.length() - 2;
    public static List<Clique> savageHaplo;
    private static AtomicInteger y = new AtomicInteger();

    private class CorrelationContainer {
        public int o11;
        public int o12;
        public int o21;
        public int o22;
        public int reads;

        public CorrelationContainer(int o11, int o12, int o21, int o22, int reads) {
            this.o11 = o11;
            this.o12 = o12;
            this.o21 = o21;
            this.o22 = o22;
            this.reads = reads;
        }
    }

    private Map<String, CorrelationContainer> correlationMap;

    /**
     * Main method for SNV method on PacBio reads input
     *
     * @param sample Given input reads with N in the beginning or end of read to make them all equal size see in {@link DataReader#readOneLined}
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public List<SNVResultContainer> getHaplotypes(IlluminaSNVSample sample) {
        long start = System.currentTimeMillis();
        System.out.println("Start 2SNV method");
        System.out.print("Compute profile");
        double[][] profile = Utils.profile(sample, al);
        //TODO remove this hack
        double t = profile[1][617];
        profile[1][617] = profile[3][617];
        profile[3][617] = t;
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.println("Profile 617: C " + sample.reads.parallelStream().filter(s -> readHasCharAtPosition(s, 617, 'C')).count() + " T "
                + sample.reads.parallelStream().filter(s -> readHasCharAtPosition(s, 617, 'T')).count());
        System.out.print("Compute split columns");
        IlluminaSNVSample splitted = DataReader.splitColumns(sample, profile, "ACGT-");
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute SNV data structure");
        SNVStructure structure = SNVStructureBuilder.buildIllumina(splitted, sample, profile);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Compute cliques");
        String consensus = Utils.consensus(structure.profile, al);
        //TODO think about it
        StringBuilder tmp = new StringBuilder(consensus);
        tmp.setCharAt(617, 'C');
        consensus = tmp.toString();
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
        System.out.println("Answers:");
        System.out.println(savageHaplo);
        //getMergedCluques first time to get cliques
        Set<Set<Integer>> cliques = run(splitted, structure, sample, true, false);
        System.out.println("Found cliques:" + cliques.stream().map(s -> new Clique(s, structure)).collect(Collectors.toList()));
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));
        System.out.print("Start getting haplotypes");
        // divide by clusters and find haplotypes
        List<SNVResultContainer> snvResultContainers = processCliques(cliques, structure, sample, false);
        System.out.println(" - DONE " + (System.currentTimeMillis() - start));

        return snvResultContainers;
    }


    public Set<Set<Integer>> run(IlluminaSNVSample splittedSample, SNVStructure struct, IlluminaSNVSample src, boolean log, boolean onlySnps) {
        correlationMap = new HashMap<>();
        int more = 0;
        String consensus = Utils.consensus(struct.profile, al);
        //TODO remove hack
        StringBuilder str = new StringBuilder(consensus);
        str.setCharAt(617, 'C');
        consensus = str.toString();
        List<Set<Integer>> adjacencyList = new ArrayList<>();
        int[][] commonReads = new int[src.referenceLength][src.referenceLength];
        System.out.println("Common reads martix calculation");

        int cores = Runtime.getRuntime().availableProcessors();
        List<Callable<Boolean>> tasks = new ArrayList<>();

        for (int i = 0; i < src.referenceLength; i++) {
            tasks.add(new ParallelTask(i, commonReads, struct, src));
        }
        ExecutorService service = Executors.newFixedThreadPool(cores);
        try {
            List<Future<Boolean>> futures = service.invokeAll(tasks);
            service.shutdown();
            futures.forEach(future -> {
                try {
                    future.get();
                } catch (InterruptedException | ExecutionException e) {
                    System.err.println("Error! Parallel tasks were not successful on get");
                    e.printStackTrace();
                }
            });
        } catch (InterruptedException e) {
            System.err.println("Error! Parallel tasks were not successful on invoke");
            e.printStackTrace();
        }
        System.out.println();
        try {
            int c = System.in.read();
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < splittedSample.referenceLength; i++) {
            adjacencyList.add(new HashSet<>());
            int l = struct.rowMinors[i].length;
            int first = i / minorCount;
            if (l < thresholdForPositions(commonReads[first][first])) {
                continue;
            }
            int[] hits = getHits(struct, struct.rowMinors[i], splittedSample.referenceLength);
            for (int j = 0; j < hits.length; j++) {
                //skip small amount of hits
                int second = j / minorCount;
                int reads = commonReads[first][second];
                if ((first != second && hits[j] >= thresholdForPositions(reads)
                        && Math.abs(first - second) >= 5)) {
                    //get unsplitted columns, minors, o_kl

                    int allele1 = getAllele(i, struct);
                    int allele2 = getAllele(j, struct);

                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    /*
                     * false 1 means that in actual sample it has another minor or N in given position
                     */
                    int o22 = hits[j];
                    int o21 = struct.rowMinors[i].length; //all 2*
                    int o12 = struct.rowMinors[j].length; //all *2
                    // subtract 2N and false 21 from o21
                    o21 = calculateO21(o21, second, i, src, struct, consensus);
                    if (o21 == 0) {
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%d p=%.3e reads=%d",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hits[j], 0.0, reads));
                        o12 = calculateO12(o12, first, j, src, struct, consensus);

                        int o11 = getO11(struct, src, consensus, i, j, first, second, o22, o21, o12, reads);
                        correlationMap.put(getCorrelationKey(first, second, m1, m2), new CorrelationContainer(o11, o12, o21, o22, reads));
                        adjacencyList.get(i).add(j);
                        continue;
                    }
                    //subtract N2 and false 12 from o12
                    o12 = calculateO12(o12, first, j, src, struct, consensus);
                    if (o12 == 0) {
                        if (log)
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%d p=%.3e reads=%d",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hits[j], 0.0, reads));
                        int o11 = getO11(struct, src, consensus, i, j, first, second, o22, o21, o12, reads);
                        correlationMap.put(getCorrelationKey(first, second, m1, m2), new CorrelationContainer(o11, o12, o21, o22, reads));
                        adjacencyList.get(i).add(j);
                        continue;
                    }


                    //subtract 1N from reads
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
                            System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%d p=%.3e reads=%d",
                                    first, second, m1, m2, l, struct.rowMinors[j].length, hits[j], p, reads));
                        correlationMap.put(getCorrelationKey(first, second, m1, m2), new CorrelationContainer(o11, o12, o21, o22, reads));
                        adjacencyList.get(i).add(j);
                    } else {
                        double pvalue = Utils.binomialPvalue(o22, p, reads);
                        if (pvalue < 0.00001 / (splittedSample.referenceLength * (splittedSample.referenceLength - 1) / 2)) {
                            if (log)
                                System.out.println(String.format("%d %d %c %c m1=%d m2=%d hits=%d p=%.3e reads=%d",
                                        first, second, m1, m2, l, struct.rowMinors[j].length, hits[j], pvalue, reads));
                            correlationMap.put(getCorrelationKey(first, second, m1, m2), new CorrelationContainer(o11, o12, o21, o22, reads));
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
        if (log) System.out.println(adjacencyList.stream().mapToInt(Set::size).sum());
        if (log) System.out.println(more);
        if (onlySnps) {
            Set<Set<Integer>> result = new HashSet<>();
            Set<Integer> positions = new HashSet<>();
            for (int i = 0; i < adjacencyList.size(); i++) {
                if (adjacencyList.get(i).size() > 0) {
                    positions.add(i);
                }
            }
            result.add(positions);
            return result;
        }
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


    private Set<Set<Integer>> getMergedCliques(List<Set<Integer>> adjacencyList, SNVStructure structure) {
        int edges = 0;
        Set<Integer> p = new HashSet<>();
        for (int i = 0; i < adjacencyList.size(); i++) {
            p.add(i);
        }
        List<Set<Integer>> cliques = AlgorithmUtils.findCliques(new HashSet<>(), p, new HashSet<>(), adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toList());
        System.out.println("Cliques before merge " + cliques.size());
        cliques.stream().map(c -> new Clique(c, structure)).forEach(System.out::println);
        List<Set<Integer>> cliquesAdjacencyMatrix = new ArrayList<>();
        cliques.forEach(i -> cliquesAdjacencyMatrix.add(new HashSet<>()));
        for (int i = 0; i < cliques.size(); i++) {
            for (int j = i + 1; j < cliques.size(); j++) {
                if (mergeScore(cliques.get(i), cliques.get(j), adjacencyList) > 0) {
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

    @Override
    protected int mergeScore(Set<Integer> c1, Set<Integer> c2, List<Set<Integer>> adjacencyList) {
        HashSet<Integer> tmp = new HashSet<>(c1);
        tmp.retainAll(c2);
        int intersection = tmp.size();
        if (intersection > 1) {
            return 1;
        }
        if (intersection < 1) {
            return 0;
        }
        int matches = 0;
        for (Integer i : c1) {
            for (Integer i2 : c2) {
                if (i / minorCount == i2 / minorCount && i % minorCount != i2 % minorCount) {
                    return 0;
                }
                if (adjacencyList.get(i).contains(i2) && !tmp.contains(i) && !tmp.contains(i2)) {
                    matches++;
                }
            }
        }
        return matches;
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
    private List<SNVResultContainer> processCliques(Set<Set<Integer>> cliques, SNVStructure struct, IlluminaSNVSample src, boolean log) {
        String z = Utils.consensus(struct.profile, al);
        //TODO remove hack
        StringBuilder tmp = new StringBuilder(z);
        tmp.setCharAt(617, 'C');
        String consensus = tmp.toString();
        //add consensus clique only if we have correlated minors with less than 40% frequency
        boolean addConsensus = false;
        outer:
        for (Set<Integer> clique : cliques) {
            for (Integer i : clique) {
                for (Integer j : clique) {
                    if (i.equals(j)) {
                        continue;
                    }
                    char m1 = al.charAt(getAllele(i, struct));
                    char m2 = al.charAt(getAllele(j, struct));
                    CorrelationContainer c = correlationMap.get(getCorrelationKey(i / minorCount, j / minorCount, m1, m2));
                    if (c != null && c.o11 / (double) c.reads < 0.4) {
                        addConsensus = true;
                        break outer;
                    }
                }
            }
        }
        //TODO think about the rule
        if (true || addConsensus) {
            cliques.add(new HashSet<>());
        }
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
                for (double[] aProfile : profile) {
                    if (aProfile[i] > max) {
                        max = aProfile[i];
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

        List<SNVResultContainer> result = haplotypes.stream().filter(distinctByKey(p -> p.haplotype)).collect(Collectors.toList());
        List<String> h = result.stream().map(s -> s.haplotype).collect(Collectors.toList());
        List<Double> frequencies = new IlluminaEM().frequencies(h, src);
        for (int i = 0; i < frequencies.size(); i++) {
            result.get(i).frequency = frequencies.get(i);
        }
        return result.stream().sorted((s1, s2) -> -Double.compare(s1.frequency, s2.frequency)).collect(Collectors.toList());
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
        class Container {
            private int distance;
            private int coincidences;

            private Container(int d, int c) {
                this.coincidences = c;
                this.distance = d;
            }
        }
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
            List<Container> distancesFromCliques = new ArrayList<>();
            for (String c : allCliquesCharacters) {
                int d = 0;
                int coincidences = 0;
                boolean isConsensusClique = c.equals(consensusClique);
                for (int i = 0; i < allPositionsInCliques.size(); i++) {
                    char charAtPosition = readCharAtPosition(read, allPositionsInCliques.get(i));
                    //increase distance
                    if (charAtPosition != 'N') {
                        if (charAtPosition != c.charAt(i)) {
                            d++;
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
                }

                if (coincidences == 0) {
                    d = 1_000_000;
                }
                distancesFromCliques.add(new Container(d, coincidences));
            }
            int minDistance = distancesFromCliques.stream().mapToInt(x -> x.distance).min().orElse(1_000_000);
            if (minDistance < 1_000_000) {
                for (int i = 0; i < distancesFromCliques.size(); i++) {
                    int d = distancesFromCliques.get(i).distance;
                    int c = distancesFromCliques.get(i).coincidences;
                    //don't add if it has only wrong coincidences
                    if (d == minDistance && d != c) {
                        clusters.get(allCliquesCharacters.get(i)).add(read);
                    }
                }
            }

        }
        return clusters;
    }

    private int thresholdForPositions(int reads) {
        return 10 * (1 + reads / 30_000);
    }

    private int calculateO12(int o12, int first, int j, IlluminaSNVSample src, SNVStructure struct, String consensus) {
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
        return o12;
    }

    private int calculateO21(int o21, int second, int i, IlluminaSNVSample src, SNVStructure struct, String consensus) {
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
        return o21;
    }

    private String getCorrelationKey(int first, int second, char m1, char m2) {
        return first + m1 + "_" + second + m2;
    }

    private class ParallelTask implements Callable<Boolean> {
        private int i;
        private int[][] commonReads;
        private SNVStructure struct;
        private IlluminaSNVSample src;

        ParallelTask(int i, int[][] commonReads, SNVStructure struct, IlluminaSNVSample src) {
            this.i = i;
            this.commonReads = commonReads;
            this.struct = struct;
            this.src = src;
        }

        @Override
        public Boolean call() throws Exception {
            System.out.print("\r" + y.incrementAndGet());
            for (int j = i; j < src.referenceLength; j++) {
                if (i == j) {
                    commonReads[i][j] = struct.readsAtPosition[i].length;
                }
                if (Math.abs(i - j) < 5) {
                    continue;
                }
                commonReads[i][j] = getCommonReadsCount(struct.readsAtPosition[i], struct.readsAtPosition[j]);
                commonReads[j][i] = commonReads[i][j];
            }
            return null;
        }
    }


}

