package by.bsu.start;

import by.bsu.algorithms.BruteForce;
import by.bsu.algorithms.BruteForceHamming;
import by.bsu.algorithms.SNVMethod;
import by.bsu.algorithms.SignatureHammingMethod;
import by.bsu.algorithms.SignatureMethod;
import by.bsu.algorithms.TreeMethod;
import by.bsu.model.IntIntPair;
import by.bsu.model.KMerDict;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.PairEndRead;
import by.bsu.model.SNVStructure;
import by.bsu.model.Sample;
import by.bsu.util.CallSignature;
import by.bsu.util.FasReader;
import by.bsu.util.SequencesTreeBuilder;
import by.bsu.util.Utils;
import by.bsu.util.builders.KMerDictBuilder;
import by.bsu.util.builders.KMerDictChunksBuilder;
import by.bsu.util.builders.SNVStructureBuilder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


public class Start {
    public static final String ACGT = "ACGT-";
    private static Path folder = Paths.get("cleaned_independent_264");
    private static List<Sample> allFiles = new ArrayList<>();
    public static Map<String, String> settings = new HashMap<>();
    public static List<Set<Integer>> snps = new ArrayList<>();

    private static void loadAllFiles(Path folder) throws IOException {
        if (allFiles.isEmpty()) {
            allFiles = FasReader.readSampleList(folder.toString(), true);
        }
    }

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        String[] strings = FasReader.readList(Paths.get("2snv/ClonesRef.fa"));
        Sample sample = FasReader.readOneLined(new File("2snv/db1/reads.fas"));
        String consensus = Utils.consensus(sample.sequences, "ACGT-N");

        for (String string : strings) {
            snps.add(new HashSet<>());
            for (int i = 0; i < string.length(); i++) {
                if (string.charAt(i) != consensus.charAt(i)){
                    snps.get(snps.size()-1).add(i);
                }
            }
        }
        String key = null;
        for (String arg : args) {
            if (arg.startsWith("-") && arg.length() > 1) {
                key = arg;
                if (key.equals("-help")) {
                    settings.put(key, "true");
                }
            } else {
                if (key == null) {
                    helpOutput(arg, true);
                }
                settings.put(key, arg);
                key = null;
            }
        }
        if (settings.get("-help") != null) {
            helpOutput(null, false);
        }
        System.out.println("Waiting for start ('q' to exit)");
        int c = System.in.read();
        if (c == 113) {
            System.exit(1);
        }

        int k = Integer.parseInt(settings.getOrDefault("-k", "10"));
        int l = Integer.parseInt(settings.getOrDefault("-l", "11"));
        folder = Paths.get(settings.getOrDefault("-dir", "cleaned_independent_264"));
        if (settings.get("-m") != null) {
            switch (settings.get("-m")) {
                case "bigData":
                    testBigDataSet(k, l);
                    break;
                case "signatureTest":
                    testDitichletAlgorithm(folder, k, l);
                    break;
                case "largeRelated":
                    testLargeRelatedSamples(folder);
                    break;
                case "snv":
                    testSNV();
                    break;
                case "assembly-snv":
                    assambling2SNV();
                    break;
                default:
                    helpOutput(settings.get("-m"), true);
            }
        } else {
            testBigDataSet(k, l);
        }
    }

    private static void assambling2SNV() {
        File dir = new File("2snv/cdc/RT_pol_SAMfiles_for_Quasirecomb");
        List<SAMRecord> reads = new ArrayList<>();
        byte re = 0;
        Map<String , List<SAMRecord>> readsSet = new HashMap<>();
        for (File file : dir.listFiles()) {
            if(!file.getName().startsWith("99_S2_L001_R1_001_2 (paired)-2 trimmed (paired) (Reads, Locally Realigned)")){
                continue;
            }
                        SamReader open = SamReaderFactory.make().open(file);
            for (SAMRecord anOpen : open) {
                reads.add(anOpen);
                if (!readsSet.containsKey(anOpen.getReadName())){
                    readsSet.put(anOpen.getReadName(), new ArrayList<>());
                }
                readsSet.get(anOpen.getReadName()).add(anOpen);
            }
        }
        List<PairEndRead> pairedReads = new ArrayList<>();
        readsSet.forEach((key, value) -> {
            if (value.size() == 1) {
                pairedReads.add(new PairEndRead(Utils.byteArrayToString(value.get(0).getReadBases()),
                        null,
                        value.get(0).getAlignmentStart(),
                        -1, key)
                );
            } else {
                SAMRecord l = value.get(0);
                SAMRecord r = value.get(0);
                pairedReads.add(new PairEndRead(Utils.byteArrayToString(l.getReadBases()),
                        Utils.byteArrayToString(l.getReadBases()),
                        l.getAlignmentStart(),
                        r.getAlignmentStart(),
                        key)
                );
            }
        });
        System.out.println();
    }

    private static void helpOutput(String arg, boolean error) {
        if (error) {
            System.out.println("Error! Wrong argument value: " + arg + " . No argument name in for of -key");
        }
        System.out.println("How to set arguments:");
        System.out.println("-k 10 -- threshold for sequence similarity(10 is default value)");
        System.out.println("-l 11 -- l-mer length for Signature method(11 is default value). Used only for equal chunks, entropy-based is used by default");
        System.out.println("-dir /usr/name/tmp/ -- folder with input. (cleaned_independent_264 is default value, except of bigData)");
        System.out.println("-outDir /usr/name/tmp/ -- folder with output.");
        System.out.println("-m bigData -- run one of predefined methods. Methods are: bigData, signatureTest. bigData is default");
        System.out.println("-algsToRun signature,tree -- which methods run for bigData method. Methods are: signature, signature-hamming, tree, brute. signature is default");
        System.out.println("-testsToRun 1,2-4,6 -- which tests to run for bigData test. Run all tests by default. Can be any combination with commas and dashes");
        System.out.println("-testsPrefix db -- which prefix do you use for all eligible tests. Then, program will run tests from folders/files db1,db2,db3,... By default program reads from folder");
        System.out.println("-threads 4 -- how many threads to use in parallel. Usually just the number of cores is the best choice");
        System.out.println("Final command can look as follows:");
        System.out.println("java -jar sequence-comparison.jar -k 10 -testsToRun 1,3-5,8 -algsToRun signature-hamming -outDir output");
        System.exit(1);
    }


    private static void testDitichletAlgorithm(Path folder, int k, int l) throws IOException, InterruptedException, ExecutionException {
        System.out.println("Start testDitichletAlgorithm with k=" + k + " l=" + l);
        long start = System.currentTimeMillis();
        loadAllFiles(folder);
        //allFiles.add(FasReader.readSampleFromFolder(new File("test_data/db1")));
        //allFiles.add(FasReader.readSampleFromFolder(new File("test_data/db2")));
        //allFiles.add(FasReader.readSampleFromFolder(new File("test_data/db3")));
        KMerDict[] kdicts = new KMerDict[allFiles.size()];
        int cores = Runtime.getRuntime().availableProcessors();
        String threads = Start.settings.get("-threads");
        if (threads != null) {
            cores = Integer.valueOf(threads);
        }
        System.out.println("Running threads = " + cores);
        ExecutorService executor = Executors.newFixedThreadPool(cores);
        List<Callable<Void>> kmerTaskList = new ArrayList<>();
        for (int j = 0; j < allFiles.size(); j++) {
            int finalJ = j;
            kmerTaskList.add(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    k1[index] = KMerDictBuilder.getDict(s1, l);
                    return null;
                }

                private Sample s1;
                private KMerDict[] k1;
                private int index;

                {
                    s1 = allFiles.get(finalJ);
                    k1 = kdicts;
                    index = finalJ;
                }
            });
        }
        List<Future<Void>> future = executor.invokeAll(kmerTaskList);
        for (Future<Void> f : future) {
            f.get();
        }
        System.out.println("Time to build dictionaries = " + (System.currentTimeMillis() - start));
        executor.shutdown();
        executor = Executors.newFixedThreadPool(cores);
        List<Callable<Long>> taskList = new ArrayList<>();
        for (int j = 0; j < allFiles.size(); j++) {
            for (int fIndex = j + 1; fIndex < allFiles.size(); fIndex++) {
                taskList.add(new CallSignature(allFiles.get(j), allFiles.get(fIndex), kdicts[j], kdicts[fIndex], k));
            }
        }
        List<Future<Long>> futures = executor.invokeAll(taskList);

        for (Future<Long> fut : futures) {
            fut.get();
        }
        executor.shutdown();
        System.out.println("testDitichletAlgorithm has ended with time = " + (System.currentTimeMillis() - start));
        System.out.println("coincidenceFilter=" + SignatureMethod.coincidenceFilter);
        System.out.println("executionCount=" + SignatureMethod.executionCount);
        System.out.println();
    }

    private static void testBigDataSet(int k, int l) throws IOException, ExecutionException, InterruptedException {
        File dir = new File(settings.getOrDefault("-dir", "test_data"));
        String[] algsToRun = settings.getOrDefault("-algsToRun", "signature").split(",");
        for (File file : dir.listFiles()) {
            if (isTestToRun(file)) {
                Sample sample = FasReader.readSampleFromFolder(file, false);
                if (Arrays.stream(algsToRun).filter(s -> s.equals("signature")).count() > 0) {
                    runSigWithTime(k, l, sample);
                }
                if (Arrays.stream(algsToRun).filter(s -> s.equals("signature-hamming")).count() > 0) {
                    runSigHamWithTime(k, l, sample);
                }
                if (Arrays.stream(algsToRun).filter(s -> s.equals("tree")).count() > 0) {
                    runTreeWithTime(k, l, sample);
                }
                if (Arrays.stream(algsToRun).filter(s -> s.equals("brute")).count() > 0) {
                    runBruteWithTime(k, sample);
                }
                if (Arrays.stream(algsToRun).filter(s -> s.equals("brute-hamming")).count() > 0) {
                    runBruteHammingWithTime(k, sample);
                }
            }
        }
    }

    private static void testSNV() throws IOException, ExecutionException, InterruptedException {
        File dir = new File(settings.getOrDefault("-dir", "2snv/db1"));
        for (File file : dir.listFiles()) {
            Sample sample = FasReader.readOneLined(file);
            long start;
            start = System.currentTimeMillis();
            double[][] profile = Utils.profile(sample, "ACGT-N");
            System.out.println("Profile time " + (System.currentTimeMillis() - start));
            Sample splitted = FasReader.splitColumns(sample, profile, "ACGT-");
            System.out.println("splitColumns " + (System.currentTimeMillis() - start));
            double[][] rotProfile = Utils.profile(splitted, "12");
            SNVStructure structure = SNVStructureBuilder.build(splitted, sample, profile);
            System.out.println("structure " + (System.currentTimeMillis() - start));
            new SNVMethod().run(splitted, structure, sample);
           System.out.println("run " + (System.currentTimeMillis() - start));
        }
    }

    private static void runSigWithTime(int k, int l, Sample query) throws ExecutionException, InterruptedException, IOException {
        long start;
        start = System.currentTimeMillis();
        KMerDict k1 = KMerDictBuilder.getDict(query, l);
        System.out.println("Dict time " + (System.currentTimeMillis() - start));
        new SignatureMethod().run(query, k1, k);
        System.out.println("Signature time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runSigHamWithTime(int k, int l, Sample query) throws ExecutionException, InterruptedException, IOException {

        long start;
        start = System.currentTimeMillis();
        double[][] profile = Utils.profile(query);
        System.out.println("Profile time " + (System.currentTimeMillis() - start));
        KMerDictChunks dict = false ? KMerDictChunksBuilder.getDict(query, l) : KMerDictChunksBuilder.getDict(query, k + 7, profile);
        new SignatureHammingMethod().runParallel(query, dict, k);
        System.out.println("Signature hamming time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runSigWithTime(int k, int l, Sample s1, Sample s2) throws ExecutionException, InterruptedException, IOException {
        long start;
        start = System.currentTimeMillis();
        KMerDict k1 = KMerDictBuilder.getDict(s1, l);
        KMerDict k2 = KMerDictBuilder.getDict(s2, l);
        long length = new SignatureMethod().run(s1, s2, k1, k2, k);
        System.out.println("count " + length);
        System.out.println("Signature time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runBruteWithTime(int k, Sample query) throws IOException {
        long start = System.currentTimeMillis();
        BruteForce.run(query, k);
        System.out.println("Brute force time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runBruteHammingWithTime(int k, Sample query) throws IOException {
        long start = System.currentTimeMillis();
        BruteForceHamming.run(query, k);
        System.out.println("Brute Hamming force time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runTreeWithTime(int k, int l, Sample query) throws IOException {
        long start = System.currentTimeMillis();
        TreeMethod.runV2(query, SequencesTreeBuilder.build(query, l), k);
        System.out.println("Tree time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static boolean isTestToRun(File file) {
        String testsToRun = settings.getOrDefault("-testsToRun", "0-1000");
        String testsPrefix = settings.getOrDefault("-testsPrefix", "db");
        List<IntIntPair> ranges = new ArrayList<>();
        for (String s : testsToRun.split(",")) {
            if (!s.contains("-")) {
                ranges.add(new IntIntPair(Integer.valueOf(s), Integer.valueOf(s)));
            } else {
                ranges.add(new IntIntPair(Integer.valueOf(s.split("-")[0]), Integer.valueOf(s.split("-")[1])));
            }
        }

        if (!file.getName().startsWith(testsPrefix)) {
            return false;
        }
        Integer testNumber = Integer.valueOf(file.getName().substring(testsPrefix.length()));
        return ranges.stream().filter(r -> r.l <= testNumber && r.r >= testNumber).count() > 0;
    }

    private static void testLargeRelatedSamples(Path folder) throws IOException {
        System.out.println("Start testLargeRelatedSamples");
        loadAllFiles(folder);
        Sample sample = allFiles.get(165);
        Sample s1 = new Sample();
        Sample s2 = new Sample();
        s1.name = "left";
        s2.name = "right";
        s1.sequences = Arrays.copyOfRange(sample.sequences, 0, sample.sequences.length / 2);
        s2.sequences = Arrays.copyOfRange(sample.sequences, sample.sequences.length / 2, sample.sequences.length);
        s1.forHamming = Utils.stringsForHamming(s1.sequences);
        s2.forHamming = Utils.stringsForHamming(s2.sequences);
        KMerDict k1 = KMerDictBuilder.getDict(s1, 11);
        KMerDict k2 = KMerDictBuilder.getDict(s2, 11);
        long start = System.currentTimeMillis();
        //System.out.println(s1.sequences.size()*s2.sequences.size());
        new SignatureMethod().run(s1, s2, k1, k2, 3);
        System.out.println("SignatureMethod time:");
        System.out.println(System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        //BruteForce.run(s1,s2, 3);
        System.out.println("BruteForce time:");
        System.out.println(System.currentTimeMillis() - start);
    }

    public static Path getOutputFilename(Sample sample, String algName) throws IOException {
        String dir = settings.getOrDefault("-outDir", "");
        return preparePath((dir.equals("") ? "" : dir + "/") + sample.name + "-" + algName + "-output.txt");
    }

    public static Path getOutputFilename(Sample sample1, Sample sample2, String algName) throws IOException {
        String dir = settings.getOrDefault("-outDir", "");
        return preparePath((dir.equals("") ? "" : dir + "/") + sample1.name + "_" + sample2.name + "-" + algName + "-output.txt");
    }

    private static Path preparePath(String filePath) throws IOException {
        Path path = Paths.get(filePath);
        Files.deleteIfExists(path);
        if (path.toFile().getParentFile() != null) {
            path.toFile().getParentFile().mkdirs();
        }
        Files.createFile(path);
        return path;
    }
}
