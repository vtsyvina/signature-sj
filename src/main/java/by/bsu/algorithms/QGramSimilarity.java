package by.bsu.algorithms;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Created by c5239200 on 4/29/17.
 * Class to calculate string q-gram similarity according to strings profile or string themself
 */
public class QGramSimilarity{

    private static final int DEFAULT_K = 3;

    private final int k;

    private static final Pattern SPACE_REG = Pattern.compile("\\s+");
    
    public QGramSimilarity(){
        k = DEFAULT_K;
    }
    
    public QGramSimilarity(int k){
        this.k = k;
    }
    
    
    
    public int similarity(Map<String, Integer> profile1, Map<String, Integer> profile2){
        int[] common = {0};
        profile1.entrySet().forEach(e-> {
            if (profile2.containsKey(e.getKey())){
                common[0]+=Math.min(e.getValue(), profile2.get(e.getKey()));
            }
        });
        return common[0];
    }
    
    public int similarity(String s1, String s2){
        return similarity(getProfile(s1), getProfile(s2));
    }

    public final Map<String, Integer> getProfile(final String string) {
        HashMap<String, Integer> shingles = new HashMap<>();

        String string_no_space = SPACE_REG.matcher(string).replaceAll(" ");
        for (int i = 0; i < (string_no_space.length() - k + 1); i++) {
            String shingle = string_no_space.substring(i, i + k);
            Integer old = shingles.get(shingle);
            if (old!=null) {
                shingles.put(shingle, old + 1);
            } else {
                shingles.put(shingle, 1);
            }
        }

        return Collections.unmodifiableMap(shingles);
    }
}
