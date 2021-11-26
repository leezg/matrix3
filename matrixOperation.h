//
// Created by lee on 2021/10/26.
//

#ifndef MATRIX_MATRIXOPERATION_H
#define MATRIX_MATRIXOPERATION_H

#include <math.h>
#include <vector>
#include <iostream>
#include <cstdio>
#define E pow(10, -12)
#define matSize 10
#define maxIterTimes 1000

using namespace std;

struct ComplexNumber {
    double Re;
    double Im;
};

double sgn(double);

const static string answerPath = "C:\\Users\\lee\\Desktop\\math\\matrix3\\answer.txt";

#endif //MATRIX_MATRIXOPERATION_H
