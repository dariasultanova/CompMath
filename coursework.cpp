#include "forsythe.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
#define N 4

double mu, rho, nu, k0, k1, k2, k3 = 24.75, x1[11];

//функция для нахождения ро
double rho_(double x) {
	return ((x * exp(x)) / pow((1 + x), 2));
}

//нелинейное уравнение для нахождения ню
double nu_(double x) {
	return (2 - x * exp(x));
}

//нелинейное уравнение для нахождения k0
double k0_(double y) {
	return (exp(y) - 2 * pow((y - 1), 2));
}

//вычисление ро, ню, k0, k1, k2, k3
void coef() {
	double errest, flag;
	int nofun;
	rho = 1.9490979 * Quanc8(rho_, 0, 1, 0.0, 1.0e-10, &errest, &nofun, &flag);

	nu = 113.7691520 * Zeroin(nu_, 0, 1, 1.0e-10);
	k0 = 225.0327352 * Zeroin(k0_, 0, 1, 1.0e-10);

	double A[3][3] = { 3.2, -1.5, 0.5,
					   1.6, 2.5, -1.0,
					   1.0, 4.1, -1.5 };

	double b[3] = { 0.9, 1.55, 2.08 };
	int ipvt[3];
	double cond = 0;

	Decomp(3, *A, &cond, ipvt);
	Solve(3, *A, b, ipvt);

	k1 = 192 * b[0];
	k2 = 19.23076923 * b[1];
	//k3 = 10 * b[2];
	return;
}

//f(x)
double f(double x) {
	if (x >= 0.2) {
		return k1 + k0;
	} else {
		return k1 + k0 + 1600;
	}
}

//P(x)
double P(double x) {
	if (x >= 0.2) {
		return (k1 + k0) * x;
	} else {
		return (k1 + k0) * x + 1600 * (x + 0.2);
	}
}

//система дифференциальных уравнений
void fun(double t, double z[N], double dz[N]) {
	dz[0] = z[2];
	dz[1] = z[3];
	dz[2] = (-f(z[0] - z[1])*(z[2] - z[3]) - mu * k3*z[2] - P(z[0] - z[1]) 
		      - k3*z[0] + k0*rho*sin(nu*t) + mu*k0*rho*nu*cos(nu*t)) / 5;
	dz[3] = (f(z[0] - z[1])*(z[2] - z[3]) - mu * k2*z[3] + P(z[0] - z[1]) 
		      - k2*z[1] - k0*rho*sin(nu*t) - mu*k0*rho*nu*cos(nu*t)) / 5;
	return;
}

//rkf45 и среднеквадратичный критерий
double F(double mu_) {
	mu = mu_;
	double sum = 0;
	double x[11] = { 10, 20.6, -18, 4.86, 3.22, -4, 1.63, 0.35, -0.81, 0.42, 0.009 };
	unsigned char work[6 * N * sizeof(Float) + sizeof(struct rkf_inside)];
	int flag = 1;
	double z[N], t = 0, re = 1.0e-5, ae = 1.0e-5, tout = 1, max = 0;
	z[0] = 10;
	z[1] = 5;
	z[2] = 100;
	z[3] = 100;
	x1[0] = z[0];
	rkf p;
	for (int i = 1; i < 11; i++) {
		p = { fun, z, t, tout, re, ae, N, flag, work };
		rkf45(&p);
		printf("t = %.0f   x1 = %4.6f\n", tout, z[0]);
		t = tout;
		tout++;
		x1[i] = z[0];
		sum += pow((x[i] - x1[i]), 2);
	}
	cout << endl;
	return sum;
}

int main() {
	coef();
	printf("%f  %f  %f  %f  %f  %f\n", rho, nu, k0, k1, k2, k3);

	mu = FMin(F, 0.01, 0.5, 1.0e-10);
	cout << "mu = " << mu << endl;

	system("pause");
}