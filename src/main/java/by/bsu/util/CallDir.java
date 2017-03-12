package by.bsu.util;

import by.bsu.algorithms.DirichletMethod;
import by.bsu.model.KMerDict;
import by.bsu.model.IntIntPair;
import by.bsu.model.Sample;

import java.util.Set;
import java.util.concurrent.Callable;

/**
 * Created by c5239200 on 2/6/17.
 */
public class CallDir implements Callable<Set<IntIntPair>> {

    private Sample sample1;
    private Sample sample2;
    KMerDict dict1;
    KMerDict dict2;
    int k;
    public CallDir(Sample sample1, Sample sample2, KMerDict dict1, KMerDict dict2, int k){
        this.sample1 = sample1;
        this.sample2 = sample2;
        this.dict1 = dict1;
        this.dict2 = dict2;
        this.k = k;
    }
    @Override
    public Set<IntIntPair> call() throws Exception {
        return DirichletMethod.run(sample1, sample2, dict1, dict2, k);
    }
}
