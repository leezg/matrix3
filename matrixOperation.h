//
// Created by lee on 2021/10/26.
//

#ifndef MATRIX_MATRIXOPERATION_H
#define MATRIX_MATRIXOPERATION_H

#include <math.h>
#include <vector>
#include <iostream>
#include <cstdio>
#include <iomanip>
#define E pow(10, -12)
#define matSize 4
#define maxIterTimes 1000
#define SIGMA pow(10, -7)

using namespace std;

struct ComplexNumber {
    double Re;
    double Im;
};

double sgn(double);
void interpolate(vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, int, int);
vector<vector<double>> fitSurface(vector<vector<double>>, int&);
void initMat(vector<vector<double>>&, int, int);

const static string answerPath = "C:\\Users\\lee\\Desktop\\math\\matrix3\\answer.txt";

#endif //MATRIX_MATRIXOPERATION_H
