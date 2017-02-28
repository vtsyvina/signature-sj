import by.bsu.algorithms.BruteForce;
import by.bsu.algorithms.DirichletMethod;
import by.bsu.algorithms.PointsMethod;
import by.bsu.model.*;
import by.bsu.util.CallDir;
import by.bsu.util.FasReader;
import by.bsu.util.KMerDictBuilder;
import by.bsu.util.PointsBuilder;
import org.apache.commons.text.beta.similarity.LevenshteinDistance;
import org.openjdk.jmh.Main;
import org.openjdk.jmh.runner.RunnerException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;


public class Start {
    static String str1 = "CACCTACAACGGCCCTGGTGGTATCGCAGTTACTCCGGATCCCACAAGCTGTCGTGGATATGGTGGCGGGGGCCCACTGGGGAGTCCTGGCGGGCCTTGCCTACTATTCCATGGTGGGGAACTGGGCTAAGGTCTTGATTGTGATGCTACTCTTTGCCGGCGTTGACGGGAGCACCCGCATGGTGGGGGAGACGGCGGGCAGGGACACCCGTGCGCTGACCAGCCTCTTCGCCCCGGGGGCGTCCCAGAAAATCCAGCTGATAC";
    static String str2 = "CCCCTACGACGGCATTGGTGGTAGCTCAGCTGCTCCGGATCCCACAAGCCATCATGGACATGATCGCTGGTGCCCACTGGGGAGTCCTGGCGGGCATAGCGTATTTCTCCATGGTGGGGAACTGGGCGAAGGTCCTGGTAGTGCTGCTGCTATTTGCCGGCGTCGACGCGGGGACCCGCGTCACCGGGGGAAGTGCCGCCTTCACCACGGCTGGGGTTGCTGGTATCTTCAGCCCAGGCGCCAAGCAGAACATCCAACTGGTCA";
    static LevenshteinDistance defaultLev = LevenshteinDistance.getDefaultInstance();
    static LevenshteinDistance limitedLev = new LevenshteinDistance(10);
    static Path folder = Paths.get("cleaned_independent_264");
    static List<Sample> allFiles = new ArrayList<>();

    static {
        for (File file : folder.toFile().listFiles()) {
            try {
                allFiles.add(new Sample(file.getName(), FasReader.readList(file.toPath())));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public static void main(String[] args) throws IOException, RunnerException, InterruptedException, ExecutionException {
//        Main.main(args);
//        testPointsAlgorithm();
        //testLargeRelatedSamples();
        //testDitichletAlgorithm();
        int cores = 4;
        List<List<Integer>> parts = new ArrayList<>();
        List<Integer> whole = new ArrayList<>();
        for (int i = 0; i < 38; i++) {
            whole.add(i);
        }
        for (int i = 0; i < cores; i++) {
            int max = i == cores-1 ? whole.size() :(i+1)*((whole.size()+1)/cores)  ;
            parts.add(whole.subList(i*((whole.size()+1)/cores),max));
        }
        testBigDataSet();
    }

    private static void testLargeRelatedSamples() {
        Sample sample = allFiles.get(165);
        Sample s1 = new Sample();
        Sample s2 = new Sample();
        s1.name = "left";
        s2.name = "right";
        s1.sequences = new HashMap<>();
        s2.sequences = new HashMap<>();
        int i =0;
        int half = sample.sequences.size()/2;
        for (Map.Entry<Integer, String> entry : sample.sequences.entrySet()){
            if (i < half){
                s1.sequences.put(entry.getKey(), entry.getValue());
            }else {
                s2.sequences.put(entry.getKey(), entry.getValue());
            }
            i++;
        }
        KMerDict k1 = KMerDictBuilder.getDict(s1, 11);
        KMerDict k2 = KMerDictBuilder.getDict(s2, 11);
        long start = System.currentTimeMillis();
        System.out.println(s1.sequences.size()*s2.sequences.size());
        List<Pair> res = DirichletMethod.run(s1, s2, k1 , k2, 3);
        System.out.println(System.currentTimeMillis()-start);
        start = System.currentTimeMillis();
        res = PointsMethod.run(s1,s2, PointsBuilder.buildPoints(s1), PointsBuilder.buildPoints(s2), 3);
        System.out.println(System.currentTimeMillis()-start);
        start = System.currentTimeMillis();
        res = BruteForce.run(s1,s2, 3);
        System.out.println(System.currentTimeMillis()-start);
    }

    private static void testPointsAlgorithm() {
        long start = System.currentTimeMillis();
        List<Points> points = new ArrayList<>();
        for (Sample sample : allFiles){
            points.add(PointsBuilder.buildPoints(sample));
        }
        long runStart1 = System.currentTimeMillis();
//        System.out.println("Start "+allFiles.get(0).name+" "+allFiles.get(40).name
//                +" pairs "+allFiles.get(0).sequences.size()*allFiles.get(40).sequences.size());
//        PointsMethod.run(allFiles.get(0), allFiles.get(40), points.get(0), points.get(40), 8);
//        System.out.println(System.currentTimeMillis() - runStart1);
        for (int j = 0; j < allFiles.size(); j++) {
            for (int k = j + 1; k < allFiles.size(); k++) {
                long runStart = System.currentTimeMillis();
                System.out.println("Start "+allFiles.get(j).name+" "+allFiles.get(k).name
                        +" pairs "+allFiles.get(j).sequences.size()*allFiles.get(k).sequences.size());
                PointsMethod.run(allFiles.get(j), allFiles.get(k), points.get(j), points.get(k), 8);
                System.out.println(System.currentTimeMillis() - runStart);
            }
        }
        System.out.println(System.currentTimeMillis() - start);
    }


    private static void testDitichletAlgorithm() throws IOException, InterruptedException, ExecutionException {

        long start = System.currentTimeMillis();
        allFiles = new ArrayList<>();
        for (File file : folder.toFile().listFiles()) {
            try {
                allFiles.add(new Sample(file.getName(), FasReader.readList(file.toPath())));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        KMerDict[] kdicts = new KMerDict[allFiles.size()];
        ExecutorService executor = Executors.newFixedThreadPool(4);
        List<Callable<Void>> kmerTaskList = new ArrayList<>();
        for (int j = 0; j < allFiles.size(); j++) {
            int finalJ = j;
            kmerTaskList.add(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    k1[index] = KMerDictBuilder.getDict(s1, 11);
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
        List<Future<Void>> future =  executor.invokeAll(kmerTaskList);
        for (Future<Void> f : future){
           f.get();
        }
        System.out.println(System.currentTimeMillis() - start);
        executor.shutdown();
        ExecutorService executor1 = Executors.newFixedThreadPool(4);
        List<Callable<List<Pair>>> taskList = new ArrayList<>();
        for (int j = 0; j < allFiles.size(); j++) {
            for (int k = j + 1; k < allFiles.size(); k++) {
                taskList.add(new CallDir(allFiles.get(j), allFiles.get(k), kdicts[j], kdicts[k], 10));
            }
        }
        List<Future<List<Pair>>> l = executor1.invokeAll(taskList);

        for (Future<List<Pair>> fut : l) {
            fut.get();
        }
        executor1.shutdown();
        System.out.println("END !!!" + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void testBigDataSet() throws IOException {
        int k = 10;
        int l = 11;
        Sample query = new Sample("query_close", FasReader.readList(Paths.get("test_data/query1/close.fas")));
        runBruteWithTime(k, query);
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);

        System.out.println();
        query = new Sample("query_medium", FasReader.readList(Paths.get("test_data/query2/medium.fas")));
        runBruteWithTime(k, query);
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);

        System.out.println();
        query = new Sample("query_medium", FasReader.readList(Paths.get("test_data/query3/far.fas")));
        runBruteWithTime(k, query);
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);

        System.out.println();
        query = new Sample("db1", FasReader.readList(Paths.get("test_data/db1/1000.fas")));
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);
        runBruteWithTime(k, query);

        System.out.println();
        query = new Sample("db2", FasReader.readList(Paths.get("test_data/db2/2000.fas")));
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);
        runBruteWithTime(k, query);

        System.out.println();
        query = new Sample("db3", FasReader.readList(Paths.get("test_data/db3/4000.fas")));
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);
        runBruteWithTime(k, query);


        System.out.println();
        query = new Sample("db4", FasReader.readList(Paths.get("test_data/db4/8000.fas")));
        runDirWithTime(k, l, query);
        runPointsWithTime(k, query);
        runBruteWithTime(k, query);
    }

    private static void runDirWithTime(int k, int l, Sample query) {
        long start;
        start = System.currentTimeMillis();
        DirichletMethod.run(query, KMerDictBuilder.getDict(query, l),k);
        System.out.println("Diri "+(System.currentTimeMillis()-start));
    }

    private static void runPointsWithTime(int k, Sample query) {
        long start;
        start = System.currentTimeMillis();
        PointsMethod.run(query, PointsBuilder.buildPoints(query), k);
        System.out.println("Points "+(System.currentTimeMillis()-start));
    }

    private static void runBruteWithTime(int k, Sample query) {
        long start = System.currentTimeMillis();
        BruteForce.run(query, k);
        System.out.println("Brute "+(System.currentTimeMillis()-start));
    }

}
