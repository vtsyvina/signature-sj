package by.bsu.algorithms;

import by.bsu.model.KMerDict;
import by.bsu.model.Sample;
import by.bsu.model.Pair;
import by.bsu.util.KMerDictBuilder;
import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.cursors.IntCursor;
import com.carrotsearch.hppc.cursors.IntIntCursor;
import com.carrotsearch.hppc.cursors.LongCursor;
import org.apache.commons.text.beta.similarity.HammingDistance;
import org.apache.commons.text.beta.similarity.LevenshteinDistance;

import java.util.*;

/**
 * Algorithm use Dirichlet method to filter pairs on sequences from sample/samples
 */
public class DirichletMethod {

        public static boolean DEBUG = true;

        public static List<Pair> run(Sample sample1, Sample sample2, KMerDict dict1, KMerDict dict2, int k){
            List<Pair> result = new ArrayList<>();
            int comps = 0;
            Set<Pair> pairsToCompare = new HashSet<>();

            int kMerCoincidences = 0;
            for (LongSet positionHashes : dict1.wholeSampleFixedPositionHashesList){
                for (LongCursor hash : positionHashes){
                    if (dict2.allHashesSet.contains(hash.value)){
                        kMerCoincidences++;
                        break;
                    }
                }
            }
            if (kMerCoincidences < dict1.fixedkMersCount - k){
                return result;
            }
            int oldLength = sample1.sequences.size()*sample2.sequences.size();
            sample1 = filterUnlikelySequences(k, dict1, dict2, sample1);
            sample2 = filterUnlikelySequences(k, dict2, dict1, sample2);
            if (sample1.sequences.size() * sample2.sequences.size() == 0){
                return result;
            }
            if (sample1.sequences.size() < sample2.sequences.size()){
                Sample tmp = sample1;
                sample1 = sample2;
                sample2 = tmp;
            }
            if (sample1.sequences.size()*sample2.sequences.size() < oldLength) {
                dict1 = KMerDictBuilder.getDict(sample1, dict1.l);
                dict2 = KMerDictBuilder.getDict(sample2, dict2.l);
            }
            boolean sameLength = sample1.sequences.values().iterator().next().length() ==
                    sample2.sequences.values().iterator().next().length();
            for ( Map.Entry<Integer, String> seqEntity : sample1.sequences.entrySet()){
                IntIntMap possibleSequences = new IntIntHashMap(dict2.sequencesNumber);
                int seq = seqEntity.getKey();
                List<Pair> tuplesToSort = new ArrayList<>();
                for (int i = 0; i < dict1.fixedkMersCount; i++) {
                    long hash = dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())[i];
                    if (dict2.allHashesSet.contains(hash)){
                        tuplesToSort.add(new Pair(i, dict2.hashToSequencesMap.get(hash).size()));
                    }
                }
                tuplesToSort.sort(Comparator.comparing(o -> o.r));
                for (int i = 0; i < tuplesToSort.size(); i++) {
                    if (i <= tuplesToSort.size() - (dict1.fixedkMersCount-k)){
                        for (IntCursor possibleSeq :
                                dict2.hashToSequencesMap
                                        .get(dict1.sequenceFixedPositionHashesList.get(seq)[tuplesToSort.get(i).l])
                                ){
                            possibleSequences.putOrAdd(possibleSeq.value, 1, 1);
                        }
                    } else {
                        IntIntMap tmp = new IntIntHashMap(possibleSequences.size());
                        for (IntIntCursor entry : possibleSequences){
                            long hash = dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())[tuplesToSort.get(i).l];
                            boolean isInSecondDIct = dict2.hashToSequencesMap.get(hash).contains(entry.key);
                            if (isInSecondDIct ||
                                    dict1.fixedkMersCount - k <= entry.value + tuplesToSort.size() - i){
                                int add = isInSecondDIct ? 1 : 0;
                                tmp.put(entry.key, entry.value + add);
                            }
                        }
                        possibleSequences = tmp;
                    }
                }
                for (IntIntCursor s : possibleSequences){
                    if (!seqEntity.getKey().equals(s.key)
                            && !pairsToCompare.contains(new Pair(seqEntity.getKey(), s.key))
                            && s.value >= dict1.fixedkMersCount - k){
                        pairsToCompare.add(new Pair(seqEntity.getKey(), s.key));
                    }
                }
            }
            LevenshteinDistance distance = new LevenshteinDistance(k);
            HammingDistance hammingDistance = new HammingDistance();
            int reduce = 0;
            for (Pair pair : pairsToCompare){
                comps++;
                if (sameLength && hammingDistance.apply(sample1.sequences.get(pair.l), sample2.sequences.get(pair.r))<=k){
                    result.add(pair);
                    reduce++;
                    continue;
                }
                int d = distance.apply(sample1.sequences.get(pair.l), sample2.sequences.get(pair.r));
                if (d != -1){
                    result.add(pair);
                }
            }
            if (DEBUG){
                System.out.println("comps = "+comps);
                System.out.println("reduce = "+reduce);
                System.out.println("length = "+result.size());
            }
            if (!result.isEmpty()){
                System.out.printf("Found %s %s%n", sample1.name, sample2.name);
            }
            return result;
        }

        private static Sample filterUnlikelySequences(int k, KMerDict dict1, KMerDict dict2, Sample sample){
            Sample result = new Sample();
            result.name = sample.name;
            result.sequences = new HashMap<>();
            for ( Map.Entry<Integer, String> seqEntity : sample.sequences.entrySet()){
                int absenceCount = 0;
                for( long hash : dict1.sequenceFixedPositionHashesList.get(seqEntity.getKey())){
                    if ( !dict2.hashToSequencesMap.containsKey(hash)){
                        absenceCount++;
                    }
                }
                if (absenceCount <= k) {
                    result.sequences.put(seqEntity.getKey(), seqEntity.getValue());
                }
            }
            return  result;
        }
}
