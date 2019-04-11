#define _USE_MATH_DEFINES
#include "Forsythe.hpp"
#include <cmath>
#include <iostream>

using namespace std;

#define h 0.2
#define ndim 11
double x;

Float fun(double t) 
{
	if (t != 0) { 
		return exp(x * t) * sin(t);
	} else { 
		return 0; 
	} 
}


double w_k(size_t k, double x, double x_[], size_t m)
{
	double res = 1;

	for (size_t i = 0; i < k; i++) {
		res *= (x - x_[i]);
	}

	for (size_t i = (k + 1); i < m; i++) {
		res *= (x - x_[i]);
	}

	return res;
}

double lagrange(double x, double x_[], size_t m, double y[])
{
	double res = 0;

	for (size_t k = 0; k < m; k++) {
		res += (y[k] * w_k(k, x, x_, m) / w_k(k, x_[k], x_, m));
	}

	return res;
}


int main() 
{
	double a = 0.0, b = 1.0;
	double ae = 0.0, re = 1.0e-10;
	double errest, flag;
	int nofun;

	double x_k[ndim], X[ndim], Y[ndim], B[ndim], C[ndim], D[ndim], S[ndim], L[ndim];

	cout << "QUANC8\n";
	cout << " x      f(x)\n";
	x = 0.0;
	for (int j = 0; j < ndim; j++) {
		X[j] = x;
		Y[j] = Quanc8(fun, a, b, ae, re, &errest, &nofun, &flag);
		x += 0.2;
		printf("%1.1f   %.6f\n", X[j], Y[j]);
	}
	cout << endl;
	

	Spline(ndim, X, Y, B, C, D);
	for (size_t k = 1; k < ndim; k++) {
		x_k[k] = (k - 0.5) * h;
		S[k] = SEval(ndim, x_k[k], X, Y, B, C, D);
	}

	
	for (size_t i = 1; i < ndim ; i++) {
		L[i] = lagrange(x_k[i], X, ndim, Y);
	}
	

	double y[ndim];
	for (size_t k = 1; k < ndim; k++) {
		y[k] = (exp(x_k[k]) * (x_k[k] * sin(1) - cos(1)) + 1) / (x_k[k] * x_k[k] + 1);
	}

	
	printf(" x      f(X)       S(X)       L(X)        (f-S)          (f-L)\n");
	for (int k = 1; k < ndim; k++) {
		printf("%1.1f   %.6f   %.6f   %.6f   %e   %e\n", 
			x_k[k], y[k], S[k], L[k], abs(y[k] - S[k]), abs(y[k] - L[k]));
	}
	system("pause");
	return 0;
}