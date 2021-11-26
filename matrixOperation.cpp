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



//曲面拟合
double *fitSurface(vector<double> z, int &kvalue) {
    int xs = 11, ys = 21, num = 9;
    vector<double> B = vector<double>(xs * num);
    vector<double> G = vector<double>(ys * num);
    vector<double> P = vector<double>(xs * ys);
    vector<double> B_temp = vector<double>(xs * num);
    vector<double> G_temp = vector<double>(ys * num);
    vector<double> B_T = vector<double>(xs * num);
    vector<double> G_T = vector<double>(ys * num);
    vector<double> BB = vector<double>(xs * num);
    vector<double> GG = vector<double>(ys * num);
    vector<double> C = vector<double>(xs * ys);

    for(int i = 0; i < num; i++) {
        for(int j = 0; j < xs; j++) {
            B[j * num + i] = pow(0.08 * j, i);
        }
        for(int j = 0; j < ys; j++) {
            G[j * num + i] = pow(0.5 + 0.05 * j, i);
        }
    }

    double sigma = 0;
    cout << endl << "选择过程的k和sigma值" << endl;
    for(int i = 0; i < num; i++) {
        sigma = 0;
        copyMat(B, xs, i+1, num, B_temp);
        transposeMat(B_temp, xs, i+1, B_T);
        multiplyMat(B_T, B_temp, i+1, xs, i+1, BB);
        inverseMat(BB, i+1, B_temp);
        copyMat(G, ys, i+1, num, G_temp);
        transposeMat(G_temp, ys, i+1, G_T);
        multiplyMat(G_T, G_temp, i+1, ys, i+1, GG);
        inverseMat(GG, i+1, G_T);
        multiplyMat(B_temp, B_T, i+1, i+1, xs, BB);
        multiplyMat(BB, z, i+1, xs, ys, GG);
        multiplyMat(GG, G_temp, i+1, ys, i+1, BB);
        multiplyMat(BB, G_T, i+1, i+1, i+1, C);
        for(int j = 0; j < xs; j++) {
            for(int k = 0; k < ys; k++) {
                double temp = 0;
                for(int p = 0; p < i+1; p++)
                    for(int q = 0; q < i+1; q++)
                        temp += C[p*(i+1)+q]*B[j*num+p]*G[k*num+q];
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