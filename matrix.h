//
// Created by lee on 2021/10/26.
//

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include "matrixOperation.h"

class Matrix {
private:
    FILE* fp;
    //存储矩阵A及R、Q、RQ矩阵
    vector<vector<double>> matrixA;
    //打印指定矩阵
    void printMatrix(vector<vector<double>>);
    //将矩阵中小于E值的数字置0
    void zeroMatrix(vector<vector<double>>&);
    //矩阵乘法
    void matrixMult(vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&);
    //初始化R、Q、RQ矩阵
    void initQR();
    //初始化矩阵A
    void initMatrixA(vector<vector<double>>&);

    void getFValue(vector<double>&, vector<double>, vector<double>);

    void getFDerivative(vector<vector<double>>&, vector<double>);

    void solveEquations(vector<vector<double>>&, vector<double>&, vector<double>&);

    double getNorm(vector<double>);

public:
    //矩阵类构造函数
    Matrix();
    //打印矩阵A
    void printMatrix();

    void Newton(vector<double>&, vector<double>)
};






#endif //MATRIX_MATRIX_H
