package by.bsu.model;

import java.util.Objects;

/**
 * Class to store (int, String) pairs
 */
public class IntStrPair {
    public int l;
    public String r;

    public IntStrPair(int l, String r) {
        this.l = l;
        this.r = r;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        IntStrPair that = (IntStrPair) o;
        return l == that.l &&
                Objects.equals(r, that.r);
    }

    @Override
    public int hashCode() {
        return Objects.hash(l, r);
    }

    @Override
    public String toString() {
        return "IntStrPair{"+ l +
                ", " + r + '}';
    }
}
