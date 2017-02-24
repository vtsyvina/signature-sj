package by.bsu.model;

import java.util.Objects;

/**
 * Container to store point coordinates
 */
public class Point {
    public int A;
    public int C;
    public int G;
    public int T;
    public int[] arr;

    public Point(int a, int c, int g, int t) {
        A = a;
        C = c;
        G = g;
        T = t;
        arr = new int[]{A,C,T,G};
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Point point = (Point) o;
        return A == point.A &&
                C == point.C &&
                G == point.G &&
                T == point.T;
    }

    @Override
    public int hashCode() {
        return Objects.hash(A, C, G, T);
    }

    @Override
    public String toString() {
        return "Point{" +
                "A=" + A +
                ", C=" + C +
                ", G=" + G +
                ", T=" + T +
                '}';
    }
}
