package com.second;

import java.util.ArrayList;

public class Test {
    private static double[][] xRes = { { -19d / 11, -12d / 11, -36d / 11 }, { -1, -1, -1 },
            { -106d / 43, -88d / 43, -78d / 43 }, { -332d / 193, -238d / 193, -216d / 193 },
            { -2d / 7, -16d / 21, -6d / 7 } };

    private static double[][] test0_a = { { 0, 2, 3 }, { 1, 2, 4 }, { 4, 5, 6 } };
    private static double[] test0_b = { 12, 17, 32 };

    private static double[][] test1_a = { { 4, 1, 1 }, { 1, 6, 1 }, { 1, 1, 8 } };
    private static double[] test1_b = { 6, 8, 10 };

    private static double[][] test2_a = { { -4, 1, 1 }, { 1, -6, 1 }, { 1, 1, -8 } };
    private static double[] test2_b = { -6, -8, -10 };

    private static double[][] test3_a = { { -4, 5, 6 }, { 7, -6, 3 }, { 6, 7, -8 } };
    private static double[] test3_b = { 6, 8, 10 };

    private static double[][] test4_a = { { 4, 3, 3 }, { 3, 6, 3 }, { 3, 3, 8 } };
    private static double[] test4_b = { 6, 8, 10 };
    
    public static ArrayList<Pair> getTests(){
        ArrayList<Pair> tests = new ArrayList<>();
        tests.add(new Pair(test0_a, test0_b));
        tests.add(new Pair(test1_a, test1_b));
        tests.add(new Pair(test2_a, test2_b));
        tests.add(new Pair(test3_a, test3_b));
        tests.add(new Pair(test4_a, test4_b));
        return tests;
    }

    static public double[][] getRes(){
        return xRes;
    }

    static public Pair createTest5(int n, double e) {
        double[][] test5_a = new double[n][n];
        double[][] temp1 = matrix.createE(n);
        for (int i = 0; i < test5_a.length; i++) {
            for (int j = i + 1; j < test5_a.length; j++) {
                temp1[i][j] = -1;
            }
        }
        double[][] temp2 = matrix.copy(temp1);
        for (int i = 0; i < temp2.length; i++) {
            for (int j = 0; j < i; j++) {
                temp2[i][j] = 1;
            }
        }
        test5_a = matrix.sum(temp1, matrix.mult(temp2, e * 2));
        double[] test5_b = new double[n];
        for (int i = 0; i < test5_b.length - 1; i++) {
            test5_b[i] = -1;
        }
        test5_b[n - 1] = 1;
        Pair test5 = new Pair(test5_a, test5_b);

        return test5;
    }
}
