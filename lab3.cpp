#include "forsythe.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

#define N 2
#define HPRINT 0.075

void f(double t, double x[N], double dx[N])
{
	dx[0] = -14 * x[0] + 13 * x[1] + cos(1 + t);
	dx[1] = 20 * x[0] - 30 * x[1] + atan(1 + t * t);
	return;
}

void rk3(void(*f)(double x, double y[], double dydx[]),
	double t, double tout, double Zn[], double h) {
	double k1[2], k2[2], k3[2];
	double fun[2], work[2];
	double tn;
	for (double t = 0; t <= tout; t += h) {
		tn = t;
		work[0] = Zn[0];
		work[1] = Zn[1];
		f(tn, work, fun);
		k1[0] = h * fun[0];
		k1[1] = h * fun[1];
		tn = t + h / 2;
		work[0] = Zn[0] + k1[0] / 2;
		work[1] = Zn[1] + k1[1] / 2;
		f(tn, work, fun);
		k2[0] = h * fun[0];
		k2[1] = h * fun[1];
		tn = t + 3 * h / 4;
		work[0] = Zn[0] + 3 * k2[0] / 4;
		work[1] = Zn[1] + 3 * k2[1] / 4;
		f(tn, work, fun);
		k3[0] = h * fun[0];
		k3[1] = h * fun[1];
		Zn[0] = Zn[0] + (2 * k1[0] + 3 * k2[0] + 4 * k3[0]) / 9;
		Zn[1] = Zn[1] + (2 * k1[1] + 3 * k2[1] + 4 * k3[1]) / 9;
	}
}

int main()
{
	unsigned char work[6 * N * sizeof(Float) + sizeof(struct rkf_inside)];
	int flag = 1;
	double x[N], t = 0, re = 0.0001, ae = 0;
	rkf p;
	x[0] = 2;
	x[1] = 0.5;
	// Решение с помощью программы rkf45
	std::cout << "For RFK45:\n";
	for (double tout = 0; tout < 1.5; tout += HPRINT)
	{
		p = { f, x, t, tout, re, ae, N, flag, work };
		rkf45(&p);
		//printf("t = %.3f   x1 = %f   x2 = %f\n", tout, x[0], x[1]);
		std::cout << std::setw(4) << std::fixed << std::setprecision(3) << "t = " << tout
			<< std::setw(10) << std::fixed << std::setprecision(6) << "x[1] = " << std::setw(9) << x[0]
			<< std::setw(10) << std::fixed << std::setprecision(6) << "x[2] = " << std::setw(9) << x[1] << '\n';
	}
	x[0] = 2;
	x[1] = 0.5;

	// Метод Рунге-Кутты 3-й степени точности
	std::cout << "\nFor Runge-Kutta method of 3rd degree witn h = 0.075:\n";
	double h = 0.075;
	for (double tout = 0; tout < 1.5; tout += HPRINT) {
		x[0] = 2;
		x[1] = 0.5;
		rk3(f, t, tout, x, h);
		std::cout << std::setw(4) << std::fixed << std::setprecision(3) << "t = " << tout
			<< std::setw(10) << std::fixed << std::setprecision(6) << "x[1] = " << std::setw(14) << x[0]
			<< std::setw(10) << std::fixed << std::setprecision(6) << "x[2] = " << std::setw(14) << x[1] << '\n';
	}

	std::cout << "\nFor Runge-Kutta method of 3rd degree witn h = 0.0075:\n";
	h = 0.0075;
	for (double tout = 0; tout < 1.5; tout += HPRINT) {
		x[0] = 2;
		x[1] = 0.5;
		rk3(f, t, tout, x, h);
		std::cout << std::setw(4) << std::fixed << std::setprecision(3) << "t = " << tout
			<< std::setw(10) << std::fixed << std::setprecision(6) << "x[1] = " << std::setw(9) << x[0]
			<< std::setw(10) << std::fixed << std::setprecision(6) << "x[2] = " << std::setw(9) << x[1] << '\n';
	}

	system("pause");
}
