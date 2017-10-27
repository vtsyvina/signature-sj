package by.bsu.util;

import com.carrotsearch.hppc.LongSet;
import com.carrotsearch.hppc.cursors.LongCursor;

import by.bsu.model.AbstractKMerDict;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

/**
 * Created by c5239200 on 7/13/17.
 */
public class AlgorithmUtils {

    /**
     * Calculates how many l-mers from fixed positions from first dictionary are in second dictionary
     */
    public static int calculateCoincidences(AbstractKMerDict dict1, AbstractKMerDict dict2) {
        int kMerCoincidences = 0;
        for (LongSet positionHashes : dict1.wholeSampleChunksHashesList) {
            for (LongCursor hash : positionHashes) {
                if (dict2.allHashesSet.contains(hash.value)) {
                    kMerCoincidences++;
                    break;
                }
            }
        }
        return kMerCoincidences;
    }

    /**
     * Using array instead of collection provide performance improvement due to low cost of iteration and random access by index
     * @param fill
     */
    public static int[] fillPossibleSequencesFromArray(int[] fill) {
        int[] tmp = new int[fill.length];
        int last = 0;
        for (int j = 0; j < fill.length; j++) {
            if (fill[j] != 0) {
                tmp[last++] = j;
            }
        }
        int[] result = new int[last];
        System.arraycopy(tmp, 0, result, 0, last);
        return result;
    }

    /**
     * Binary search return index of x in sorted array or index of first element that is lower than x.
     * Returns 0 if x is less then arr[0]
     * @param arr Sorted array
     * @param x Value to find
     * @return index of x
     */
    public static int binarySearch(int arr[], int x)
    {
        int l = 0, r = arr.length - 1;
        while (l <= r)
        {
            int m = l + (r-l)/2;

            // Check if x is present at mid
            if (arr[m] == x)
                return m;

            // If x greater, ignore left half
            if (arr[m] < x)
                l = m + 1;

                // If x is smaller, ignore right half
            else
                r = m - 1;
        }

        // if we reach here, then element was not present
        return r == -1 ? 0 : r;
    }

    /**
     * Simple implementation of Bron - Kerbosch algorithm for finding all maximum cliques (psedocode is on Wiki)
     */
    public static Set<Set<Integer>> findCliques(Set<Integer> clique, Set<Integer> p, Set<Integer> x, List<Set<Integer>> adjacencyList) {
        Set<Set<Integer>> result = new HashSet<>();
        if (p.isEmpty() && x.isEmpty()) {
            if (clique.size() > 1){
                result.add(clique);
            }

        }
        Iterator<Integer> it = p.iterator();
        while (it.hasNext()) {
            Integer v = it.next();
            Set<Integer> newCLique = new HashSet<>(clique);
            newCLique.add(v);
            Set<Integer> pIntersection = new HashSet<>(p);
            Set<Integer> xIntersection = new HashSet<>(x);
            pIntersection.retainAll(adjacencyList.get(v));
            xIntersection.retainAll(adjacencyList.get(v));
            result.addAll(findCliques(newCLique, pIntersection, xIntersection, adjacencyList));
            it.remove();
            x.add(v);
        }
        return result;
    }

    public static int[] connectedComponents(List<Set<Integer>> adjacencyList){
        int[] result = new int[adjacencyList.size()];
        int count = 0;
        int currentComponent = 1;
        Queue<Integer> q = new LinkedList<>();
        q.add(0);
        result[0] = 1;
        while (count < result.length){
            if (q.isEmpty()){
                for (int i = 0; i < result.length; i++) {
                    if (result[i] == 0){
                        q.add(i);
                        currentComponent++;
                        result[i] = currentComponent;
                        break;
                    }
                }
            }
            Integer poll = q.poll();
            count++;
            for (Integer v : adjacencyList.get(poll)) {
                if (result[v] == 0){
                    result[v] = result[poll];
                    q.add(v);
                }
            }
        }
        return result;
    }
}
