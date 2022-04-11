package com.second;

import java.util.ArrayList;

public class matrix {

    static private int perest = 0;

    static public double[][] sum(double[][] a, double[][] b) {
        double[][] c = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                c[i][j] = a[i][j] + b[i][j];
            }
        }
        return c;
    }

    static public double[] sum(double[] a, double[] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < c.length; i++) {
            c[i] = a[i] + b[i];
        }
        return c;
    }

    static public double[][] mult(double[][] a, double[][] b) {
        double[][] c = new double[a.length][b[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                for (int k = 0; k < b.length; k++)
                    c[i][j] += a[i][k] * b[k][j];
            }
        }
        return c;
    }

    static public double[] mult(double[][] a, double[] b) {
        double[] c = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                c[i] += a[i][j] * b[j];
            }
        }
        return c;
    }

    static public double[] mult(double[] a, double b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            c[i] = b * a[i];
        }
        return c;
    }

    static public double[][] mult(double[][] a, double b) {
        double[][] c = new double[a.length][a.length];
        for (int i = 0; i < c.length; i++) {
            for (int j = 0; j < c.length; j++) {
                c[i][j] = a[i][j] * b;
            }
        }
        return c;
    }

    static public double[] mult(double[] a, double[][] b) {
        double[] c = new double[a.length];
        for (int i = 0; i < c.length; i++) {
            double temp = 0;
            for (int j = 0; j < c.length; j++) {
                temp += a[j] * b[i][j];
            }
            c[i] = temp;
        }
        return c;
    }

    static public double mult(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < b.length; i++) {
            s += a[i] * b[i];
        }
        return s;
    }

    static public double[][] transp(double[][] a) {
        double[][] c = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                c[i][j] = a[j][i];
            }
        }
        return c;
    }

    static public double norma(double[][] a) {
        double max = 0, s = 0;
        for (int i = 0; i < a.length; i++) {
            s = 0;
            for (int j = 0; j < a[0].length; j++) {
                s += Math.abs(a[i][j]);
            }
            if (s > max)
                max = s;
        }
        return max;
    }

    static public double norma(double[] a) {
        double res = 0;
        for (int i = 0; i < a.length; i++) {
            res += a[i] * a[i];
        }
        return Math.sqrt(res);
    }

    static public double norma2(double[][] a) {
        return Math.sqrt(matrix.maxEig(matrix.mult(matrix.transp(a), a)));
    }

    static public double[][] copy(double[][] a) {
        double[][] b = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[0].length; j++)
                b[i][j] = a[i][j];
        return b;
    }

    static public double[] copy(double[] a) {
        double[] c = new double[a.length];
        for (int i = 0; i < c.length; i++) {
            c[i] = a[i];
        }
        return c;
    }

    static public double determ(double[][] a) {
        double det = 1;
        matrix.perest = 0;
        double[][] c = new double[a.length][a.length];
        c = matrix.copy(a);
        c = matrix.uptringle(c);
        for (int i = 0; i < a.length; i++) {
            double prod = c[i][i];
            det *= prod;
        }
        if (matrix.perest % 2 == 1)
            det *= -1;
        matrix.perest = 0;
        return det;
    }

    static public double[][] uptringle(double[][] a) {
        int n = a.length;
        int f = a[0].length;

        double[][] c = new double[n][f];

        c = matrix.copy(a);

        for (int k = 0; k < n; k++) {

            if (c[k][k] == 0) {
                int l = 0;
                for (int h = k; h < n; h++)
                    if (c[h][k] != 0) {
                        l = h;
                        break;
                    }
                if (c[l][k] != 0) {
                    matrix.perest += 1;
                    c = matrix.perestanovka(c, k, l);
                } else
                    continue;
            }

            double del = c[k][k];

            for (int i = k + 1; i < n; i++) {
                double prod = c[i][k];
                for (int j = k; j < f; j++) {

                    c[i][j] -= c[k][j] * prod / del;
                }
            }
        }
        return c;
    }

    static public double[][] createE(int n) {
        double[][] E = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                E[i][j] = 0f;

                if (i == j)
                    E[i][j] = 1f;
            }
        return E;
    }

    static public double[][] inverse(double[][] a) {
        double[][] c = matrix.copy(a);
        int n = c.length;

        double[][] E = matrix.createE(n);

        for (int k = 0; k < n; k++) {

            if (c[k][k] == 0) {
                int l = 0;
                for (int h = k; h < n; h++)
                    if (c[h][k] != 0) {
                        l = h;
                        break;
                    }
                if (c[l][k] != 0) {
                    c = matrix.perestanovka(c, k, l);
                    E = matrix.perestanovka(E, k, l);
                } else
                    continue;
            }

            double del = c[k][k];
            for (int i = 0; i < n; i++) {
                c[k][i] /= del;
                E[k][i] /= del;
            }

            for (int i = k + 1; i < n; i++) {
                double prod = c[i][k];
                for (int j = 0; j < n; j++) {

                    c[i][j] -= c[k][j] * prod;
                    E[i][j] -= E[k][j] * prod;
                }
            }
        }

        for (int k = n - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                double prod = c[i][k];

                for (int j = 0; j < n; j++) {
                    c[i][j] -= c[k][j] * prod;
                    E[i][j] -= E[k][j] * prod;
                }
            }
        }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = E[i][j];

        return c;
    }

    static public double[][] perestanovka(double[][] a, int k, int l) {
        double[][] c = matrix.copy(a);
        for (int m = 0; m < a[0].length; m++) {
            double o, p;
            p = c[k][m];
            o = c[l][m];
            c[k][m] = o;
            c[l][m] = p;
        }
        return c;
    }

    static public void printMatrix(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                System.out.print(a[i][j] + " ");
            }
            System.out.println();
        }
    }

    static public void printVector(double[] a) {
        for (int i = 0; i < a.length; i++) {
            System.out.println(a[i]);
        }
    }

    static public int mainEl(double[][] a, int n) {
        int max = n;
        for (int i = n + 1; i < a.length; i++) {
            if (Math.abs(a[max][n]) < Math.abs(a[i][n]))
                max = i;
        }

        return max;
    }

    static public boolean isMainDiag(double[][] a) {
        double[] delta = new double[a.length];
        for (int i = 0; i < delta.length; i++) {
            double temp = 0;
            for (int j = 0; j < a.length; j++)
                temp += (i != j) ? Math.abs(a[i][j]) : 0;
            delta[i] = Math.abs(a[i][i]) - temp;
        }

        double min = delta[0];
        for (int i = 1; i < delta.length; i++) {
            if (min > delta[i])
                min = delta[i];
        }
        return (min > 0) ? true : false;
    }

    static public boolean isPositDef(double[][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a.length; j++) {
                if (a[i][j] != a[j][i])
                    return false;
            }
        }

        for (int i = 1; i < a.length; i++) {
            double[][] temp = new double[i][i];
            for (int j = 0; j < temp.length; j++)
                for (int k = 0; k < temp.length; k++)
                    temp[j][k] = a[j][k];
            if (matrix.determ(temp) < 0)
                return false;
        }
        return true;
    }

    static public double maxEig(double[][] a) {
        int n = a.length;
        double e = 10e-6;

        double[][] A = matrix.copy(a);

        double[] xk = new double[n];
        double[] xk1 = new double[n];
        xk1[0] = 1;

        double l1 = 1, l2 = 0;

        while (Math.abs(l1 - l2) > e) {
            xk = xk1;
            xk1 = matrix.mult(A, xk);
            l1 = l2;
            l2 = matrix.mult(xk1, xk) / matrix.mult(xk, xk);
        }

        return l2;
    }

    static public ArrayList<double[][]> QandR(double[][] a) {
        int n = a.length;

        ArrayList<double[][]> Q = new ArrayList<>();
        Q.add(matrix.createE(n));

        ArrayList<double[][]> QR = new ArrayList<>();

        double[][] R = matrix.copy(a);

        for (int k = 0; k < n - 1; k++) {

            double[] y = new double[n - k];
            for (int i = 0; i < y.length; i++) {
                y[i] = R[i + k][k];
            }
            double alpha = (-1) * matrix.norma(y);
            double[] z = new double[n - k];
            z[0] = 1;

            double r = matrix.norma(matrix.sum(y, matrix.mult(z, alpha)));
            double[] w = matrix.mult(matrix.sum(y, matrix.mult(z, alpha)), 1 / r);

            double[][] wwt = new double[n - k][n - k];
            for (int i = 0; i < n - k; i++) {
                for (int j = 0; j < n - k; j++) {
                    wwt[i][j] = w[i] * w[j];
                }
            }

            double[][] Qtemp = new double[n - k][n - k];
            double[][] E = matrix.createE(n - k);
            double[][] Qn = matrix.createE(n);
            double[][] R1 = new double[n - k][n - k];

            Qtemp = matrix.sum(E, matrix.mult(wwt, -2));

            for (int i = 0; i < R1.length; i++) {
                for (int j = 0; j < R1.length; j++) {
                    R1[i][j] = R[i + k][j + k];
                }
            }

            R1 = matrix.mult(Qtemp, R1);

            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    Qn[i][j] = Qtemp[i - k][j - k];
                    R[i][j] = R1[i - k][j - k];
                }
            }

            Q.add(Qn);
        }

        double[][] Qres = matrix.createE(n);
        for (int i = 0; i < Qres.length; i++) {
            Qres = matrix.mult(Qres, Q.get(i));
        }

        QR.add(Qres);
        QR.add(R);

        return QR;
    }

}
