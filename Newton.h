//
// Created by lee on 2021/11/28.
//

#ifndef MATRIX3_NEWTON_H
#define MATRIX3_NEWTON_H

#include "matrixOperation.h"

class Newton {
private:
    void getFValue(vector<double>&, vector<double>, vector<double>);

    void getFDerivative(vector<vector<double>>&, vector<double>);

    void solveEquations(vector<vector<double>>&, vector<double>&, vector<double>&);

    double getNorm(vector<double>);

public:
    //矩阵类构造函数
    Newton();

    void newton(vector<double>&, vector<double>);
};


#endif //MATRIX3_NEWTON_H
