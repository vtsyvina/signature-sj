package by.bsu.start;

import by.bsu.algorithms.BruteForce;
import by.bsu.algorithms.BruteForceHamming;
import by.bsu.algorithms.SignatureHammingMethod;
import by.bsu.algorithms.SignatureMethod;
import by.bsu.model.KMerDict;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;
import by.bsu.util.DataReader;
import by.bsu.util.Utils;
import by.bsu.util.builders.KMerDictBuilder;
import by.bsu.util.builders.KMerDictChunksBuilder;
import by.bsu.util.tasks.CallEditSignature;
import by.bsu.util.tasks.CallHammingSignature;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


public class Start {
    public static Map<String, String> settings = new HashMap<>();


    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
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

        int k = Integer.parseInt(settings.getOrDefault("-k", "10"));
        int l = Integer.parseInt(settings.getOrDefault("-l", "11"));
        File input = new File(settings.getOrDefault("-in", "cleaned_independent_264/AMC_P01_1b.fas"));
        if (settings.get("-m") != null) {
            switch (settings.get("-m")) {
                case "edit-single":
                    runEditDistance(input, k, l);
                    break;
                case "hamming-single":
                    runHammingDistance(input, k, l);
                    break;
                case "edit-multi":
                    runMulti(input, k, l, true);
                    break;
                case "hamming-multi":
                    runMulti(input, k, l, false);
                    break;
                case "brute-edit-single":
                    runBruteWithTime(k, input);
                    break;
                default:
                    helpOutput(settings.get("-m"), true);
            }
        } else {
            helpOutput(null, false);
        }
    }

    private static void helpOutput(String arg, boolean error) {
        //TODO rewrite this
        if (error) {
            System.out.println("Error! Wrong argument value: " + arg + " . No argument name in for of -key");
        }
        System.out.println("How to set arguments:");
        System.out.println("-k 10 -- threshold for sequence similarity(10 is default value)");
        System.out.println("-l 11 -- l-mer length for Signature method(11 is default value). Used only for equal chunks, entropy-based is used by default");
        System.out.println("-dir /usr/name/tmp/ -- folder with input. (cleaned_independent_264 is default value, except of bigData)");
        System.out.println("-in /usr/name/tmp/reads.fas -- input file(2snv/realigned/reads.fas is default for snv method)");
        System.out.println("-outDir /usr/name/tmp/ -- folder with output.");
        System.out.println("-m bigData -- getMergedCluques one of predefined methods. Methods are: bigData, signatureTest. bigData is default");
        System.out.println("-algsToRun signature,tree -- which methods getMergedCluques for bigData method. Methods are: signature, signature-hamming, tree, brute. signature is default");
        System.out.println("-testsToRun 1,2-4,6 -- which tests to getMergedCluques for bigData test. Run all tests by default. Can be any combination with commas and dashes");
        System.out.println("-testsPrefix db -- which prefix do you use for all eligible tests. Then, program will getMergedCluques tests from folders/files db1,db2,db3,... By default program reads from folder");
        System.out.println("-threads 4 -- how many threads to use in parallel. Usually just the number of cores is the best choice");
        System.out.println("Final command can look as follows:");
        System.out.println("java -jar sequence-comparison.jar -k 10 -testsToRun 1,3-5,8 -algsToRun signature-hamming -outDir output");
        System.exit(1);
    }

    private static void runEditDistance(File file, int k, int l) throws IOException {
        Sample sample = getSample(file);
        if (sample == null) return;
        long start = System.currentTimeMillis();
        KMerDict dict = KMerDictBuilder.getDict(sample, l);
        new SignatureMethod().runParallel(sample, dict, k);
        System.out.println("Total run time: " + (System.currentTimeMillis() - start) + ", ms");
    }


    private static void runHammingDistance(File file, int k, int l) throws IOException {
        Sample sample = getSample(file);
        long start = System.currentTimeMillis();
        double[][] profile = Utils.profile(sample);
        KMerDictChunks dict = KMerDictChunksBuilder.getDict(sample, k + 7, profile);
        new SignatureHammingMethod().runParallel(sample, dict, k);
        System.out.println("Total run time: " + (System.currentTimeMillis() - start) + ", ms");
    }


    private static void runMulti(File folder, int k, int l, boolean edit) throws IOException, InterruptedException, ExecutionException {
        if (!folder.isDirectory()){
            System.out.println("Input is not a directory");
            return;
        }
        List<Sample> samples = DataReader.readSampleList(folder, true);
        String method = edit ? "edit distance" : "hamming distance";
        System.out.println("Start "+method+" multi run for t=" + k);
        System.out.println("Total number of sample = "+samples.size());
        long start = System.currentTimeMillis();
        KMerDict[] merDicts = new KMerDict[samples.size()];
        KMerDictChunks[] chunkDicts = new KMerDictChunks[samples.size()];
        int cores = Runtime.getRuntime().availableProcessors();
        String threads = Start.settings.get("-threads");
        if (threads != null) {
            cores = Integer.valueOf(threads);
        }
        System.out.println("Running threads = " + cores);
        ExecutorService executor = Executors.newFixedThreadPool(cores);
        List<Callable<Void>> kmerTaskList = new ArrayList<>();
        for (int j = 0; j < samples.size(); j++) {
            int finalJ = j;
            kmerTaskList.add(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    if (edit){
                        mers[index] = KMerDictBuilder.getDict(s1, l);
                    } else {
                        //TODO think about how to run with entropy-based chunks
                        chunks[index] = KMerDictChunksBuilder.getDict(s1, l);
                    }
                    return null;
                }

                private Sample s1;
                private KMerDict[] mers;
                private int index;
                private KMerDictChunks[] chunks;

                {
                    s1 = samples.get(finalJ);
                    mers = merDicts;
                    chunks = chunkDicts;
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
        for (int j = 0; j < samples.size(); j++) {
            for (int fIndex = j + 1; fIndex < samples.size(); fIndex++) {
                if (edit){
                    taskList.add(new CallEditSignature(samples.get(j), samples.get(fIndex), merDicts[j], merDicts[fIndex], k));
                } else {
                    taskList.add(new CallHammingSignature(samples.get(j), samples.get(fIndex), chunkDicts[j], chunkDicts[fIndex], k));
                }

            }
        }
        List<Future<Long>> futures = executor.invokeAll(taskList);

        for (Future<Long> fut : futures) {
            fut.get();
        }
        executor.shutdown();
        System.out.println(method+" run has ended with time = " + (System.currentTimeMillis() - start));
        System.out.println();
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
        KMerDictChunks dict = KMerDictChunksBuilder.getDict(query, k + 7, profile);
        new SignatureHammingMethod().runParallel(query, dict, k);
        System.out.println("Signature hamming time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runBruteWithTime(int k, File file) throws IOException {
        Sample sample = getSample(file);
        long start = System.currentTimeMillis();
        BruteForce.run(sample, k);
        System.out.println("Brute force time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void runBruteHammingWithTime(int k, Sample query) throws IOException {
        long start = System.currentTimeMillis();
        BruteForceHamming.run(query, k);
        System.out.println("Brute Hamming force time " + (System.currentTimeMillis() - start));
        System.out.println();
    }

    public static Path getOutputFilename(Sample sample, String algName) throws IOException {
        String dir = settings.getOrDefault("-outDir", "output/");
        return preparePath((dir.equals("") ? "" : dir + "/") + sample.name + "-" + algName + "-output.txt");
    }

    public static Path getOutputFilename(Sample sample1, Sample sample2, String algName) throws IOException {
        String dir = settings.getOrDefault("-outDir", "output/");
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

    private static Sample getSample(File file) throws IOException {
        if (!file.exists()) {
            System.out.println(String.format("Input file %s does not exists", file.getCanonicalPath()));
            return null;
        }
        Sample sample;
        if (file.isDirectory()) {
            sample = DataReader.readSampleFromFolder(file, false);
        } else {
            sample = DataReader.readSampleFromFile(file);
        }
        return sample;
    }
}
