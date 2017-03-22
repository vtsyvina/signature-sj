package by.bsu.util;

import by.bsu.model.IntIntPair;
import by.bsu.model.Sample;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class that provides operations with Sample
 */
public class SampleTransformer {
    public static Sample reduceEqualLetters(Sample src){
        Sample result = new Sample();
        result.name = src.name+"-reducedLetters";
        result.sequences = new HashMap<>();
        String first = src.sequences.values().iterator().next();
        int l = first.length();
        List<IntIntPair> toStay = new ArrayList<>();
        int start = -1;
        for (int i = 0; i <= l; i++) {
            boolean fl = true;
            if (i == l){//to cover edge cases
                fl = false;
            } else {
                for (String seq : src.sequences.values()) {
                    if (seq.charAt(i) != first.charAt(i)) {
                        fl = false;
                        break;
                    }
                }
            }
            if (!fl && start == -1){
                start = i;
            } else if (fl && start != -1){
                toStay.add(new IntIntPair(start, i-1));
                start = -1;
            }
        }
        if (toStay.isEmpty()){
            toStay.add(new IntIntPair(0, l-1));
        }
        int sum = 0;
        for(IntIntPair pair: toStay){
            sum += pair.r-pair.l +1;
        }
        for (Map.Entry<Integer, String > entry : src.sequences.entrySet()){
            StringBuilder builder = new StringBuilder();
            for (IntIntPair pair : toStay){
                builder.append(entry.getValue().substring(pair.l, pair.r+1));
            }
            result.sequences.put(entry.getKey(), builder.toString());
        }
        return result;
    }
}
