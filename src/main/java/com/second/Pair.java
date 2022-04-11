package com.second;

public class Pair {
    private double[] first1;
    private int second1;

    private double[][] matr;
    private double[] vector;

    public Pair(double[] a, int count)
    {
        this.first1 = matrix.copy(a);
        this.second1 = count;
    }

    public Pair(double[][] a, double[] b){
        this.matr = matrix.copy(a);
        this.vector = matrix.copy(b);
    }

    public double[] getX(){
        return this.first1;
    }

    public int getCount(){
        return this.second1;
    }

    public double[][] getMatr(){
        return this.matr;
    }

    public double[] getVect(){
        return this.vector;
    }
}
