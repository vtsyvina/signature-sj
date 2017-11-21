package by.bsu.util.tasks;

import by.bsu.algorithms.SignatureHammingMethod;
import by.bsu.model.KMerDictChunks;
import by.bsu.model.Sample;

import java.util.concurrent.Callable;

public class CallHammingSignature implements Callable<Long> {
    private Sample sample1;
    private Sample sample2;
    private KMerDictChunks dict1;
    private KMerDictChunks dict2;
    private int k;

    public CallHammingSignature(Sample sample1, Sample sample2, KMerDictChunks dict1, KMerDictChunks dict2, int k) {
        this.sample1 = sample1;
        this.sample2 = sample2;
        this.dict1 = dict1;
        this.dict2 = dict2;
        this.k = k;
    }

    @Override
    public Long call() throws Exception {
        return new SignatureHammingMethod().run(sample1, sample2, dict1, dict2, k);
    }
}
