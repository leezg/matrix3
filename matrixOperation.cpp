//
// Created by lee on 2021/10/26.
//

#include "matrixOperation.h"

void initMat(vector<vector<double>>& mat, int m, int n) {
    for (int i = 0; i < m; i++) {
        mat.push_back(vector<double>(n));
    }
}

vector<vector<double>> transposeMat(vector<vector<double>> A) {
    vector<vector<double>> A_T;
    for (int i = 0; i < A[0].size(); i++) {
        A_T.push_back(vector<double>(A.size()));
    }
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            A_T[j][i] = A[i][j];
        }
    }
    return A_T;
}

void matrixMult(vector<vector<double>> ma, vector<vector<double>> mb, vector<vector<double>>& ms)  {
    ms = vector<vector<double>>();
    for (int i = 0; i < ma.size(); i++) {
        ms.push_back(vector<double>(mb[i].size()));
        for (int j = 0; j < mb[i].size(); j++) {
            ms[i][j] = 0;
            for (int k = 0; k < mb.size(); k++) {
                ms[i][j] += ma[i][k] * mb[k][j];
            }
        }
    }
}

void inverseMat(vector<vector<double>>& A, vector<vector<double>>& result)
{
    vector<vector<double>> A_temp = A;
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++) {
            if(i != j) {
                result[i][j] = 0;
            } else {
                result[i][j] = 1;
            }
        }

    for(int i = 0; i < A.size(); i++) {
        int index = i;
        for (int j = i + 1; j < A.size(); j++) {
            if (abs(A_temp[j][i]) > abs(A_temp[index][i])) {
                index = j;
            }
        }

        if (index != i) {
            for(int j = 0; j < A.size(); j++) {
                double temp = A_temp[index][j];
                A_temp[index][j] = A_temp[i][j];
                A_temp[i][j] = temp;
                temp = result[index][j];
                result[index][j] = result[i][j];
                result[i][j] = temp;
            }
        }

        double temp = A_temp[i][i];
        if (temp != 1) {
            for(int j = 0; j < A.size(); j++) {
                A_temp[i][j] /= temp;
                result[i][j] /= temp;
            }
        }
        for (int j = 0; j < A.size(); j++) {
            if (A_temp[j][i] != 0 && i != j) {
                temp = A_temp[j][i];
                for (int k = 0; k < A.size(); k++) {
                    A_temp[j][k] -= temp * A_temp[i][k];
                    result[j][k] -= temp * result[i][k];
                }
            }
        }
    }
}

//曲面拟合
vector<vector<double>> fitSurface(vector<vector<double>> z, int& kvalue) {
    int xs = 11, ys = 21, num = 9;
    vector<vector<double>> B, G, P, C;
    for (int i = 0; i < num; i++) {
        B.push_back(vector<double>(xs));
        G.push_back(vector<double>(ys));
    }
    for (int i = 0; i < xs; i++) {
        P.push_back(vector<double>(ys));
    }

    for (int i = 0; i < num; i++) {
        for (int j = 0; j < xs; j++) {
            B[i][j] = pow(0.08 * j, i);
        }
        for (int j = 0; j < ys; j++) {
            G[i][j] = pow(0.5 + 0.05 * j, i);
        }
    }

    double sigma = 0;
    cout << endl << "选择过程的k和sigma值" << endl;
    for (int i = 0; i < num; i++) {
        sigma = 0;
        vector<vector<double>> B_temp;
        vector<vector<double>> G_temp;
        for (int j = 0; j < i + 1; j++) {
            B_temp.push_back(B[j]);
            G_temp.push_back(G[j]);
        }

        vector<vector<double>> B_T = transposeMat(B_temp);
        vector<vector<double>> G_T = transposeMat(G_temp);

        vector<vector<double>> BB, GG;
        matrixMult(B_temp, B_T, BB);
        matrixMult(G_temp, G_T, GG);

        inverseMat(BB, B_T);
        inverseMat(GG, G_temp);

        vector<vector<double>> subB_T;
        for (int j = 0; j <= i; j++) {
            subB_T.push_back(B_T[j]);
        }

        matrixMult(subB_T, B_temp, BB);
        matrixMult(BB, z, GG);
        matrixMult(GG, G_T, BB);

        vector<vector<double>> subG_temp;
        for (int j = 0; j <= i; j++) {
            subG_temp.push_back(vector<double>(i + 1));
            for (int k = 0; k <= i; k++) {
                subG_temp[j][k] = G_temp[j][k];
            }
        }

        matrixMult(BB, subG_temp, C);

        int kkkk = 0;
        for(int j = 0; j < xs; j++) {
            for(int k = 0; k < ys; k++) {
                double temp = 0;
                for(int p = 0; p < i + 1; p++) {
                    for(int q = 0; q < i + 1; q++) {
                        //TODO: maybe bug
                        temp += C[p][q] * B[p][j] * G[q][k];
                    }

                }

                P[j][k] = temp;
                sigma += (z[j][k] - temp) * (z[j][k] - temp);
            }
        }
        cout << "k = " << i << setprecision(11) << setiosflags(ios::scientific|ios::uppercase) << " sigma = " << sigma << endl;
        if(sigma < SIGMA) {
            kvalue = i;
            cout << endl << "达到精度要求时的k和sigma值" << endl;
            cout << "k = " << i << setprecision(11) << setiosflags(ios::scientific|ios::uppercase) << " sigma = " << sigma << endl;
            cout << endl << "p(x, y)中的系数Crs" << endl;
            for(int k = 0; k <= i; k++) {
                for(int j = 0; j <= i; j++) {
                    cout << "C[" << k << "][" << j << "] = " << C[k][j] << endl;
                }
            }
            return C;
        }
    }
    return vector<vector<double>>();
}

void interpolate(vector<vector<double>>& t, vector<vector<double>>& u, vector<vector<double>>& z, int xs, int ys) {
    vector<double> tt = {0, 0.2, 0.4, 0.6, 0.8, 1};
    vector<double> uu = {0, 0.4, 0.8, 1.2, 1.6, 2};
    vector<vector<double>> zz = {{-0.5, -0.34, 0.14, 0.94, 2.06, 3.5},
                                 {-0.42, -0.5, -0.26, 0.3, 1.18, 2.38},
                                 {-0.18, -0.5, -0.5, -0.18, 0.46, 1.42},
                                 {0.22, -0.34, -0.58, -0.5, -0.1, 0.62},
                                 {0.78, -0.02, -0.5, -0.66, -0.5, -0.02},
                                 {1.5, 0.46, -0.26, -0.66, -0.74, -0.5}};
    double h = 0.2, tao = 0.4;
    int n = 5, m = 5;
    for(int i = 0; i < xs; i++) {
        for(int j = 0; j < ys; j++) {

            int xp = 0, yp = 0;
            if(t[i][j] <= tt[1] + h / 2) {
                xp = 1;
            }

            else if(t[i][j] > tt[n-1]-h/2) {
                xp = n-1;
            }

            else {
                for(int q = 2; q <= n-2; q++) {
                    if(t[i][j] > tt[q] - h / 2 && t[i][j] <= tt[q] + h / 2) {
                        xp = q;
                    }
                }

            }
            if(u[i][j] <= uu[1] + tao / 2) {
                yp = 1;
            } else if(u[i][j] > uu[m - 1] - tao / 2) {
                yp = n-1;
            } else {
                for(int q = 2; q <= m-2; q++) {
                    if((u[i][j] > uu[q] - tao / 2) && (u[i][j] <= uu[q] + tao / 2)) {
                        yp = q;
                    }
                }
            }
            z[i][j] = 0;
            for(int k = xp - 1; k <= xp + 1; k++) {
                double lk = 1;
                for(int ti = xp-1; ti <= xp + 1; ti++) {
                    if(ti != k) {
                        lk *= (t[i][j] - tt[ti]) / (tt[k] - tt[ti]);
                    }
                }
                for(int r = yp - 1; r <= yp + 1; r++) {
                    double lr = 1;
                    for(int ti = yp - 1; ti <= yp + 1; ti++) {
                        if(ti != r)
                            lr *= (u[i][j] - uu[ti]) / (uu[r] - uu[ti]);
                    }
                    z[i][j] += lk * lr * zz[k][r];
                }
            }
        }
    }
}
