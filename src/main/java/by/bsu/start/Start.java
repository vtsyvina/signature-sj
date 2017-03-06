package by.bsu.start;

import by.bsu.algorithms.BruteForce;
import by.bsu.algorithms.DirichletMethod;
import by.bsu.algorithms.PointsMethod;
import by.bsu.algorithms.TreeMethod;
import by.bsu.model.*;
import by.bsu.model.Pair;
import by.bsu.util.*;
import com.sun.tools.javac.util.*;
import org.openjdk.jmh.Main;
import org.openjdk.jmh.runner.RunnerException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;

import static java.util.stream.Collectors.toMap;


public class Start {
    private static Path folder = Paths.get("cleaned_independent_264");
    private static List<Sample> allFiles = new ArrayList<>();

    static {
//        for (File file : folder.toFile().listFiles()) {
//            try {
//                allFiles.add(new Sample(file.getName(), FasReader.readList(file.toPath())));
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
    }

    public static void main(String[] args) throws IOException, RunnerException, InterruptedException, ExecutionException {
//        Main.main(args);
//        testPointsAlgorithm();
        //testLargeRelatedSamples();
        //testDitichletAlgorithm();
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
        Set<Pair> res = DirichletMethod.run(s1, s2, k1 , k2, 3);
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
//        System.out.println("by.bsu.start.Start "+allFiles.get(0).name+" "+allFiles.get(40).name
//                +" pairs "+allFiles.get(0).sequences.size()*allFiles.get(40).sequences.size());
//        PointsMethod.run(allFiles.get(0), allFiles.get(40), points.get(0), points.get(40), 8);
//        System.out.println(System.currentTimeMillis() - runStart1);
        for (int j = 0; j < allFiles.size(); j++) {
            for (int k = j + 1; k < allFiles.size(); k++) {
                long runStart = System.currentTimeMillis();
                System.out.println("by.bsu.start.Start "+allFiles.get(j).name+" "+allFiles.get(k).name
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
        List<Callable<Set<Pair>>> taskList = new ArrayList<>();
        for (int j = 0; j < allFiles.size(); j++) {
            for (int k = j + 1; k < allFiles.size(); k++) {
                taskList.add(new CallDir(allFiles.get(j), allFiles.get(k), kdicts[j], kdicts[k], 10));
            }
        }
        List<Future<Set<Pair>>> l = executor1.invokeAll(taskList);

        for (Future<Set<Pair>> fut : l) {
            fut.get();
        }
        executor1.shutdown();
        System.out.println("END !!!" + (System.currentTimeMillis() - start));
        System.out.println();
    }

    private static void testBigDataSet() throws IOException, ExecutionException, InterruptedException {
        int k = 3;
        int l = 11;
        Sample query = new Sample("query_close", FasReader.readList(Paths.get("test_data/query1/close.fas")));
        runBruteWithTime(k, query);
        runTreeWithTime(k, query);
        //runDirWithTime(k, l, query);
        //runPointsWithTime(k, query);

        System.out.println();
        query = new Sample("query_medium", FasReader.readList(Paths.get("test_data/query2/medium.fas")));
        runBruteWithTime(k, query);
        //runDirWithTime(k, l, query);
        //runPointsWithTime(k, query);

        System.out.println();
        query = new Sample("query_far", FasReader.readList(Paths.get("test_data/query3/far.fas")));
        runBruteWithTime(k, query);
        //runDirWithTime(k, l, query);
        //runPointsWithTime(k, query);

        System.out.println();
        query = new Sample("db1", FasReader.readList(Paths.get("test_data/db1/1000.fas")));
        //runDirWithTime(k, l, query);
        //runTreeWithTime(k, query);

        System.out.println();
        query = new Sample("db2", FasReader.readList(Paths.get("test_data/db2/2000.fas")));
        //runDirWithTime(k, l, query);
        //runTreeWithTime(k, query);

        System.out.println();
        query = new Sample("db3", FasReader.readList(Paths.get("test_data/db3/4000.fas")));
        //runDirWithTime(k, l, query);
        //runTreeWithTime(k, query);


        System.out.println();
        query = new Sample("db4", FasReader.readList(Paths.get("test_data/db4/8000.fas")));
        runDirWithTime(k, l, query);
        //runTreeWithTime(k, query);

        System.out.println();
        query = new Sample("db5", FasReader.readList(Paths.get("test_data/db5/16000.fas")));
        //runDirWithTime(k, l, query);
        //runTreeWithTime(k, query);

        System.out.println();
        query = new Sample("db6", FasReader.readList(Paths.get("test_data/db6/32000.fas")));
        //runDirWithTime(k, l, query);

        System.out.println();
        Map<Integer, String > seq = FasReader.readList(Paths.get("test_data/db7/32000 (1).fas"));
        Map<Integer, String > tmp = FasReader.readList(Paths.get("test_data/db7/32000 (2).fas"));
        int[] size = new int[1];
        size[0] = seq.size();
        tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
        seq.putAll(tmp);
        query = new Sample("db7", seq);
        runDirWithTime(k, l, query);

        System.out.println();
        seq = FasReader.readList(Paths.get("test_data/db8/32000.fas"));
        tmp = FasReader.readList(Paths.get("test_data/db8/32000 (2).fas"));
        size[0] = seq.size();
        tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
        seq.putAll(tmp);
        tmp = FasReader.readList(Paths.get("test_data/db8/32000 (3).fas"));
        size[0] = seq.size();
        tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
        seq.putAll(tmp);
        tmp = FasReader.readList(Paths.get("test_data/db8/32000 (4).fas"));
        size[0] = seq.size();
        tmp = tmp.entrySet().stream().collect(toMap(e -> e.getKey() + size[0], Map.Entry::getValue));
        seq.putAll(tmp);
        query = new Sample("db8", seq);

        runDirWithTime(k, l, query);
    }

    private static void runDirWithTime(int k, int l, Sample query) throws ExecutionException, InterruptedException {
        long start;
        start = System.currentTimeMillis();
        KMerDict k1 = KMerDictBuilder.getDict(query, l);
        Set<Pair> r = DirichletMethod.run(query, k1 ,k);
        System.out.println("Diri "+(System.currentTimeMillis()-start));
        r.clear();
    }

    private static void runPointsWithTime(int k, Sample query) {
        long start;
        start = System.currentTimeMillis();
        PointsMethod.run(query, PointsBuilder.buildPoints(query), k);
        System.out.println("Points "+(System.currentTimeMillis()-start));
    }

    private static void runBruteWithTime(int k, Sample query) {
        long start = System.currentTimeMillis();
        Set<Pair> r = BruteForce.run(query, k);
        System.out.println("Brute "+(System.currentTimeMillis()-start));
    }

    private static void runTreeWithTime(int k, Sample query) {
        long start = System.currentTimeMillis();
        TreeMethod.run(query, SequencesTreeBuilder.build(query), k);
        System.out.println("Tree "+(System.currentTimeMillis()-start));
    }
}
