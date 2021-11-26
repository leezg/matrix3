//
// Created by lee on 2021/11/11.
//

#include "matrix.h"

void Matrix::initQR() {
    for (int i = 0; i < maxLength; i++) {
        matrixQ.push_back(vector<double>());
        matrixR.push_back(vector<double>());
        matrixQR.push_back(vector<double>());
        for (int j = 0; j < maxLength; j++) {
            if (i == j) {
                matrixQ[i].push_back(1);
            } else {
                matrixQ[i].push_back(0);
            }
            matrixR[i].push_back(matrixA[i][j]);
            matrixQR[i].push_back(0);
        }
    }
    zeroMatrix(matrixR);
}

void Matrix::getQR() {
    initQR();
    for(int r = 0; r < maxLength - 1; r++){
        double d = 0;
        double c = 0;
        double h = 0;
        vector<double> vectorU;
        vector<double> vectorW;
        vector<double> vectorP;
        for(int i = r + 1; i < maxLength; i++) {
            d += matrixR[i][r] * matrixR[i][r];
        }
        if(abs(d) == 0) {
            continue;
        }
        d += matrixR[r][r] * matrixR[r][r];
        d = sqrt(d);
        if (matrixR[r][r] == 0) {
            c = d;
        }
        else {
            c = -sgn(matrixR[r][r]) * d;
        }
        h = c * c - c * matrixR[r][r];
        for(int i = 0; i < r; i++) {
            vectorU.push_back(0);
        }
        vectorU.push_back(matrixR[r][r] - c);
        for(int i = r + 1; i < maxLength; i++) {
            vectorU.push_back(matrixR[i][r]);
        }
        for(int i = 0; i < maxLength; i++) {
            vectorW.push_back(0);
            for(int j = 0; j < maxLength; j++) {
                vectorW[i] += matrixQ[i][j] * vectorU[j];
            }
        }
        for(int i = 0; i < maxLength; i++) {
            for(int j = 0; j < maxLength; j++)
                matrixQ[i][j] -= vectorW[i] * vectorU[j] / h;
        }
        for(int i = 0; i < maxLength; i++){
            vectorP.push_back(0);
            for(int j = 0; j < maxLength; j++) {
                vectorP[i] += matrixR[j][i] * vectorU[j] / h;
            }
        }
        for(int i = 0; i < maxLength; i++) {
            for (int j = 0; j < maxLength; j++)
                matrixR[i][j] -= vectorU[i] * vectorP[j];
        }
    }
    zeroMatrix(matrixR);
    zeroMatrix(matrixQ);
    matrixMult(matrixR, matrixQ, matrixQR);
    zeroMatrix(matrixQR);
//    cout << "Q" << endl;
    fprintf(fp, "Q,\n");
    printMatrix(matrixQ);
//    cout << "R" << endl;
    fprintf(fp, "R,\n");
    printMatrix(matrixR);
//    cout << "RQ" << endl;
    fprintf(fp, "RQ\n");
    printMatrix(matrixQR);
}

void Matrix::QRMethod() {
    vector<vector<double>> matrixM;
    vector<ComplexNumber> L = vector<ComplexNumber>(10);
    double det = 0;
    for (int i = 0; i < maxLength; i++) {
        matrixM.push_back(vector<double>());
        for (int j = 0; j < maxLength; j++) {
            matrixM[i].push_back(0);
        }
    }

    int m = maxLength - 1;
    int r = 0;
    for (int k = 0; k < maxIterTimes; k++) {
        if (m == 0) {
            L[r].Re = matrixA[m][m];
            L[r].Im = 0;
            break;
        } else if (m < 0) {
            break;
        }
        if (abs(matrixA[m][m - 1]) < E) {
            L[r].Re = matrixA[m][m];
            L[r].Im = 0;
            m--;
            r++;
        } else {
            if (m == 1) {
                det = (matrixA[m][m] + matrixA[m - 1][m - 1]) * (matrixA[m][m] + matrixA[m - 1][m - 1])
                        - 4 * (matrixA[m][m] * matrixA[m - 1][m - 1] - matrixA[m - 1][m] * matrixA[m][m - 1]);
                if (det > 0) {
                    L[r].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2 + sqrt(det) / 2;
                    L[r].Im = 0;
                    L[r + 1].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2 - sqrt(det) / 2;
                    L[r + 1].Im = 0;
                } else {
                    L[r].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2;
                    L[r].Im = sqrt(-det) / 2;
                    L[r + 1].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2;
                    L[r + 1].Im = -sqrt(-det) / 2;
                }
                m -= 2;
                r += 2;
                continue;
            } else if (abs(matrixA[m - 1][m - 2]) < E) {
                det = (matrixA[m][m] + matrixA[m - 1][m - 1]) * (matrixA[m][m] + matrixA[m - 1][m - 1])
                        - 4 * (matrixA[m][m] * matrixA[m - 1][m - 1] - matrixA[m - 1][m] * matrixA[m][m - 1]);
                if (det > 0) {
                    L[r].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2 + sqrt(det) / 2;
                    L[r].Im = 0;
                    L[r + 1].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2 - sqrt(det) / 2;
                    L[r + 1].Im = 0;
                }
                else {
                    L[r].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2;
                    L[r].Im = sqrt(-det) / 2;
                    L[r + 1].Re = (matrixA[m][m] + matrixA[m - 1][m - 1]) / 2;
                    L[r + 1].Im = -sqrt(-det) / 2;
                }
                m -= 2;
                r += 2;
                continue;
            } else {
                double s;
                double t;
                s = matrixA[m - 1][m - 1] + matrixA[m][m];
                t = matrixA[m - 1][m - 1] * matrixA[m][m] - matrixA[m][m - 1] * matrixA[m - 1][m];
                matrixMult(matrixA, matrixA, matrixM);
                for (int i = 0; i < 10; i++)
                {
                    for (int j = 0; j < 10; j++) {
                        matrixM[i][j] -= s * matrixA[i][j];
                    }
                    matrixM[i][i] += t;
                }
                iterate(matrixM, m + 1);
                zeroMatrix(matrixA);
            }
        }
    }
    zeroMatrix(matrixA);
    fprintf(fp, "after QR method\n");
    printMatrix(matrixA);

    for (int r = 0; r < 10; r++) {
        fprintf(fp, "\n");
        if (L[r].Im == 0) {
            fprintf(fp, "lambda[%d] = (%.12e + i*%.12e)\n", r + 1, L[r].Re, L[r].Im);
            gauss(L[r].Re);
        }
        else {
            fprintf(fp, "lambda[%d] = (%.12e + i*%.12e)\n", r + 1, L[r].Re, L[r].Im);
        }
    }
}

void Matrix::iterate(vector<vector<double>> &matrixM, int m) {
    vector<double> vectorU = vector<double>(maxLength);
    vector<double> vectorV = vector<double>(maxLength);
    vector<double> vectorP = vector<double>(maxLength);
    vector<double> vectorQ = vector<double>(maxLength);
    vector<double> vectorW = vector<double>(maxLength);
    for (int r = 0; r < m - 1; r++) {
        zeroMatrix(matrixM);
        double d = 0;
        for (int i = r + 1; i < m; i++) {
            d += matrixM[i][r] * matrixM[i][r];
        }
        if (abs(d) == 0) {
            continue;
        } else {
            double c;
            double h;
            d += matrixM[r][r] * matrixM[r][r];
            d = sqrt(d);
            if (matrixM[r][r] != 0) {
                c = -abs(matrixM[r][r]) / matrixM[r][r] * d;
            } else {
                c = d;
            }
            h = c * c - c * matrixM[r][r];
            for (int i = 0; i < r; i++) {
                vectorU[i] = 0;
            }
            vectorU[r] = matrixM[r][r] - c;
            for (int i = r + 1; i < m; i++) {
                vectorU[i] = matrixM[i][r];
            }
            for (int i = 0; i < m; i++) {
                vectorV[i] = 0;
                for (int j = 0; j < m; j++) {
                    vectorV[i] += matrixM[j][i] * vectorU[j] / h;
                }
            }
            for (int i = 0; i < m; i++) {
                vectorP[i] = 0;
                vectorQ[i] = 0;
                for (int j = 0; j < m; j++) {
                    matrixM[i][j] -= vectorU[i] * vectorV[j];
                    vectorP[i] += matrixA[j][i] * vectorU[j] / h;
                    vectorQ[i] += matrixA[i][j] * vectorU[j] / h;
                }
            }
            double t = 0;
            for (int i = 0; i < m; i++) {
                t += vectorP[i] * vectorU[i] / h;
            }
            for (int i = 0; i < m; i++) {
                vectorW[i] = vectorQ[i] - t * vectorU[i];
            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    matrixA[i][j] -= (vectorW[i] * vectorU[j] + vectorU[i] * vectorP[j]);
                }
            }
        }
    }
}