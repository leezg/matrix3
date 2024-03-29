//
// Created by lee on 2021/11/28.
//

#include "Newton.h"

Newton::Newton() {

}

//获得F的值的向量
void Newton::getFValue(vector<double>& b, vector<double> offset, vector<double> x) {
    b[1] = -(0.5 * cos(x[1]) + x[2] + x[3] + x[4] - offset[1]);
    b[2] = -(x[1] + 0.5 * sin(x[2]) + x[3] + x[4] - offset[2]);
    b[3] = -(0.5 * x[1] + x[2] + cos(x[3]) + x[4] - offset[3]);
    b[4] = -(x[1] + 0.5 * x[2] + x[3] + sin(x[4]) - offset[4]);
}

//获得F的导数矩阵
void Newton::getFDerivative(vector<vector<double>>& Fd, vector<double> x) {
    Fd[1][1] = -0.5*sin(x[1]);
    Fd[1][2] = 1;
    Fd[1][3] = 1;
    Fd[1][4] = 1;

    Fd[2][1] = 1;
    Fd[2][2] = 0.5*cos(x[2]);
    Fd[2][3] = 1;
    Fd[2][4] = 1;

    Fd[3][1] = 0.5;
    Fd[3][2] = 1;
    Fd[3][3] = -sin(x[3]);
    Fd[3][4] = 1;

    Fd[4][1] = 1;
    Fd[4][2] = 0.5;
    Fd[4][3] = 1;
    Fd[4][4] = cos(x[4]);
}

void Newton::solveEquations(vector<vector<double>>& mat, vector<double>& b, vector<double>& x)
{
    for(int k = 1; k < matSize; k++) {
        int i = k;
        for(int j = k; j <= matSize; j++) {
            if(mat[j][k] > mat[i][k])
                i = j;
        }
        for(int j = k; j <= matSize; j++) {
            double temp = mat[k][j];
            mat[k][j] = mat[i][j];
            mat[i][j] = temp;
        }
        double temp = b[k];
        b[k] = b[i];
        b[i] = temp;
        for(int i = k+1; i <= matSize; i++) {
            double mik = mat[i][k]/mat[k][k];
            for(int j = k+1; j <= matSize; j++)
                mat[i][j] = mat[i][j]-mik*mat[k][j];
            b[i] = b[i]-mik*b[k];
        }
    }
    x[matSize] = b[matSize]/mat[matSize][matSize];
    for(int k = matSize-1; k >= 1; k--) {
        x[k] = 0;
        for(int j = k+1; j <= matSize; j++)
            x[k] -= mat[k][j]*x[j];
        x[k] += b[k];
        x[k] = x[k] / mat[k][k];
    }
}

double Newton::getNorm(vector<double> vec) {
    double temp = 0;
    for(int i = 0; i < matSize; i++) {
        if(temp < fabs(vec[i])) {
            temp = fabs(vec[i]);
        }

    }
    return temp;
}

//Newton法解非线性方程组
void  Newton::newton(vector<double>& x, vector<double> offset) {
    vector<double> delta_x = vector<double>(matSize + 1);
    vector<double> tempx = vector<double>(matSize + 1);
    vector<double> F = vector<double>(matSize + 1);
    vector<vector<double>> FDerivative;
    for (int i = 0; i < matSize + 1; i++) {
        FDerivative.push_back(vector<double>(matSize + 1));
    }

    for (int i = 1; i <= matSize; i++) {
        x[i] = 1;
    }

    for (int i = 1; i < maxIterTimes; i++) {
        getFDerivative(FDerivative, x);
        getFValue(F, offset, x);
        solveEquations(FDerivative, F, delta_x);
        if(getNorm(delta_x) / getNorm(x) <= E) {
            break;
        } else {
            for(int i = 1; i <= matSize; i++) {
                x[i] += delta_x[i];
            }
        }
    }
}