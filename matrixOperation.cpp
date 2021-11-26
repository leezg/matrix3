//
// Created by lee on 2021/10/26.
//

#include "matrixOperation.h"

double sgn(double n) {
    if (n > 0) {
        return 1;
    } else if (n < 0) {
        return -1;
    } else {
        return 0;
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
}

void matrixMult(vector<vector<double>>& ma, vector<vector<double>>& mb, vector<vector<double>>& ms)  {
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

//曲面拟合
double *fitSurface(vector<vector<double>> z, int &kvalue) {
    int xs = 11, ys = 21, num = 9;
    vector<vector<double>> B, G, P, C;
    for (int i = 0; i < num; i++) {
        B.push_back(vector<double>(xs));
        G.push_back(vector<double>(ys));
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
        for (int j = 0; j < i; j++) {
            B_temp.push_back(B[j]);
            G_temp.push_back(G[j]);
        }

        vector<vector<double>> B_T = transposeMat(B_temp);
        vector<vector<double>> G_T = transposeMat(G_temp);

        vector<vector<double>> BB, GG;
        matrixMult(B_T, B_temp, BB);
        matrixMult(G_T, G_temp, GG);

        inverseMat(BB, i+1, B_temp);
        inverseMat(GG, i+1, G_T);

        matrixMult(B_temp, B_T, BB);
        matrixMult(BB, z, GG);
        matrixMult(GG, G_temp, BB);
        matrixMult(BB, G_T, C);

        for(int j = 0; j < xs; j++) {
            for(int k = 0; k < ys; k++) {
                double temp = 0;
                for(int p = 0; p < i+1; p++) {
                    for(int q = 0; q < i+1; q++) {
                        temp += C[p * (i + 1) + q] * B[j * num + p] * G[k * num + q];
                    }

                }

                P[j*ys+k] = temp;
                sigma += (z[j*ys+k]-temp)*(z[j*ys+k]-temp);
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
                    cout << "C[" << k << "][" << j << "] = " << C[k*(i+1)+j] << endl;
                }
            }
            return C;
        }
    }
    return NULL;
}