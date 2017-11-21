package by.bsu.util.tasks;

import by.bsu.algorithms.SignatureMethod;
import by.bsu.model.KMerDict;
import by.bsu.model.Sample;

import java.util.concurrent.Callable;

/**
 * Created by c5239200 on 2/6/17.
 */
public class CallEditSignature implements Callable<Long> {

    private Sample sample1;
    private Sample sample2;
    private KMerDict dict1;
    private KMerDict dict2;
    private int k;
    public CallEditSignature(Sample sample1, Sample sample2, KMerDict dict1, KMerDict dict2, int k){
        this.sample1 = sample1;
        this.sample2 = sample2;
        this.dict1 = dict1;
        this.dict2 = dict2;
        this.k = k;
    }
    @Override
    public Long call() throws Exception {
        return new SignatureMethod().run(sample1, sample2, dict1, dict2, k);
    }
}
