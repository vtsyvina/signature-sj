package by.bsu.model;

/**
 * Class to store (int, int) pairs with values from range [0, 1 000 000)
 */
public class IntIntPair {
    public int l;
    public int r;
    private int hash = -1;

    public IntIntPair(int left, int right){
        this.l = left;
        this.r = right;
    }

    @Override
    public int hashCode() {
        if (hash != -1)
            return hash;
        hash = l *1_000_000+ r;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        IntIntPair p = (IntIntPair)obj;
        return p.l == this.l && p.r == this.r;
    }

    @Override
    public String toString() {
        return "("+ l +", "+ r +")";
    }
}
