package by.bsu.model;

/**
 * Class to represent unordered pair of two positive integers up to 1 million
 * (x, y).equals((y,x))
 * l automatically >= r (for improved performance)
 */
public class Pair {
    public int l;
    public int r;
    private int hash = -1;

    public Pair(int left, int right){
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
        Pair p = (Pair)obj;
        return p.l == this.l && p.r == this.r;
    }

    @Override
    public String toString() {
        return "("+ l +", "+ r +")";
    }
}
