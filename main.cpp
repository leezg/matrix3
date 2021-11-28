#include "Newton.h"

int main() {
    system("chcp 65001");
    vector<vector<double>> t, u, z;
    vector<double> offset = vector<double>(matSize + 1);
    vector<double> x = vector<double>(matSize + 1);
    initMat(t, 11, 21);
    initMat(u, 11, 21);
    initMat(z, 11, 21);

    Newton newton = Newton();
    int kvalue = 0;
    for (int i = 0; i <= 10; i++) {
        for (int j = 0; j <= 20; j++) {
            offset[1] = 2.67 + 0.08 * i;
            offset[2] = 1.07 + 0.5 + 0.05 * j;
            offset[3] = 3.74 + 0.08 * i;
            offset[4] = 0.79 + 0.5 + 0.05 * j;
            newton.newton(x, offset);
            t[i][j] = x[1];
            u[i][j] = x[2];
        }
    }
    interpolate(t, u, z, 11, 21);
    cout << "数表: (xi, yi, f(xi, yi))" << endl;
    for (int i = 0; i <= 10; i++) {
        for (int j = 0; j <= 20; j++) {
            cout << resetiosflags(ios::scientific | ios::uppercase) << "x = " << 0.08*i << "   y = " << 0.5+0.05*j;
            cout << setprecision(11) << setiosflags(ios::scientific | ios::uppercase) << "   f(x, y) = " << z[i][j] << endl;
        }
    }
    vector<vector<double>> C = fitSurface(z, kvalue);
    for (int i = 1; i <= 8; i++) {
        for (int j = 1; j <= 5; j++) {
            offset[1] = 2.67 + 0.1 * i;
            offset[2] = 1.07 + 0.5 + 0.2 * j;
            offset[3] = 3.74 + 0.1 * i;
            offset[4] = 0.79 + 0.5 + 0.2 * j;
            newton.newton(x, offset);
            t[i-1][j-1] = x[1];
            u[i-1][j-1] = x[2];
        }
    }
    interpolate(t, u, z, 8, 5);
    double p[8][5];
    for(int i = 1; i <= 8; i++) {
        for(int j = 1; j <= 5; j++) {
            double temp = 0;
            for(int ii = 0; ii <= kvalue; ii++) {
                for(int jj = 0; jj <= kvalue; jj++) {
                    temp += C[ii][jj] * pow(0.1 * i, ii) * pow(0.5 + 0.2 * j, jj);
                }
            }
            p[i-1][j-1] = temp;
        }
    }

    cout << endl << "数表: (xi*, yi*, f(xi*, yi*), p(xi*, yi*))" << endl;
    for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 5; j++) {
            cout << resetiosflags(ios::scientific | ios::uppercase);
            cout << "x[" << i << "] = " << 0.1*i << " y[" << j << "] = " << 0.5+0.2*j << endl;
            cout << setprecision(11) << setiosflags(ios::scientific | ios::uppercase);
            cout << "p(x, y) = " << p[i][j] << "  f(x, y) = " << z[i][j] << endl;
            cout << "delta = " << p[i][j] - z[i][j] << endl;
        }
    }
    return 0;
    return 0;
}
