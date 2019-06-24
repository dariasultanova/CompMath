#define _USE_MATH_DEFINES
#include "forsythe.hpp"
#include <cmath>
#include <iostream>

int main() {
	//n = 4
	const int dim = 4;
	double C[dim][dim], Ct[dim][dim], CtC[dim][dim];
	double d[dim], d1[dim];
	int ipvt[dim];
	double cond, delta;
	double m = 1, n = 1;

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			C[i][j] = 1 / (m + n - 1);
			n++;
		}
		n = 1; m++;
	}

	m = 1; n = 1;
	for (int i = 0; i < dim; i++) {
		d[i] = 0;
		for (int j = 0; j < dim; j++) {
			d[i] = 1 / (m + n - 1) + d[i];
			n++;
		}
		n = 1; m++;
	}

	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			Ct[i][j] = C[j][i];

	m = 1; n = 1;
	for (int i = 0; i < dim; i++) {
		d1[i] = 0;
		for (int j = 0; j < dim; j++) {
			d1[i] = d1[i] + Ct[i][j] * d[j];
			n++;
		}
		n = 1; m++;
	}

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			CtC[i][j] = 0;
			for (int k = 0; k < dim; k++)
				CtC[i][j] = CtC[i][j] + Ct[i][k] * C[k][j];
		}
	}

	printf("      Matrix C               d\n");
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++)
			printf("%1.2f  ", C[i][j]);
		printf("|  %1.2f", d[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim, *C, &cond, ipvt);
	printf("cond = %e\n\n", cond);
	Solve(dim, *C, d, ipvt);
	int z = 1;
	for (int i = 0; i < dim; i++) {
		printf("x1[%d] =  %1.2f", z, d[i]);
		z++;
		printf("\n");
	}
	printf("_______________________________\n\n");
	
	printf("      Matrix CtC           Ct*d\n");
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++)
			printf("%1.2f  ", CtC[i][j]);
		printf("|  %1.2f", d1[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim, *CtC, &cond, ipvt);
	printf("cond = %e\n\n", cond);
	Solve(dim, *CtC, d1, ipvt);
	z = 1;
	for (int i = 0; i < dim; i++) {
		printf("x2[%d] =  %1.2f", z, d1[i]);
		z++;
		printf("\n");
	}
	printf("_______________________________\n");

	double s1 = 0, s2 = 0;
	for (int i = 0; i < dim; i++) {
		s1 = pow(d[i], 2) + s1;
		s2 = pow((d[i] - d1[i]), 2) + s2;
	}
	s1 = sqrt(s1); s2 = sqrt(s2);
	delta = s2 / s1;
	printf("delta = %e\n\n", delta);
	printf("===============================================================================\n");

	// n = 6
	const int dim1 = 6;
	double C1[dim1][dim1], Ct1[dim1][dim1], CtC1[dim1][dim1];
	double d_1[dim1], d1_1[dim1];
	int ipvt1[dim1];
	double cond1, delta1;
	double m1 = 1, n1 = 1;

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim1; j++) {
			C1[i][j] = 1 / (m1 + n1 - 1);
			n1++;
		}
		n1 = 1; m1++;
	}

	m1 = 1; n1 = 1;
	for (int i = 0; i < dim1; i++) {
		d_1[i] = 0;
		for (int j = 0; j < dim1; j++) {
			d_1[i] = 1 / (m1 + n1 - 1) + d_1[i];
			n1++;
		}
		n1 = 1; m1++;
	}

	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim1; j++)
			Ct1[i][j] = C1[j][i];

	m1 = 1; n1 = 1;
	for (int i = 0; i < dim1; i++) {
		d1_1[i] = 0;
		for (int j = 0; j < dim1; j++) {
			d1_1[i] = d1_1[i] + Ct1[i][j] * d_1[j];
			n1++;
		}
		n1 = 1; m1++;
	}

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim1; j++) {
			CtC1[i][j] = 0;
			for (int k = 0; k < dim1; k++)
				CtC1[i][j] = CtC1[i][j] + Ct1[i][k] * C1[k][j];
		}
	}

	printf("      Matrix C                          d\n");
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim1; j++)
			printf("%1.2f  ", C1[i][j]);
		printf("|  %1.2f", d_1[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim1, *C1, &cond1, ipvt1);
	printf("cond = %e\n\n", cond1);
	Solve(dim1, *C1, d_1, ipvt1);
	int z1 = 1;
	for (int i = 0; i < dim1; i++) {
		printf("x1[%d] =  %1.2f", z1, d_1[i]);
		z1++;
		printf("\n");
	}
	printf("___________________________________________\n\n");

	printf("      Matrix CtC                       Ct*d\n");
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim1; j++)
			printf("%1.2f  ", CtC1[i][j]);
		printf("|  %1.2f", d1_1[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim1, *CtC1, &cond1, ipvt1);
	printf("cond = %e\n\n", cond1);
	Solve(dim1, *CtC1, d1_1, ipvt1);
	z1 = 1;
	for (int i = 0; i < dim1; i++) {
		printf("x2[%d] =  %1.2f", z1, d1_1[i]);
		z1++;
		printf("\n");
	}
	printf("_______________________________\n");

	double s1_1 = 0, s2_1 = 0;
	for (int i = 0; i < dim1; i++) {
		s1_1 = pow(d_1[i], 2) + s1_1;
		s2_1 = pow((d_1[i] - d1_1[i]), 2) + s2_1;
	}
	s1_1 = sqrt(s1_1);
	s2_1 = sqrt(s2_1);
	delta1 = s2_1 / s1_1;
	printf("delta = %e\n\n", delta1);
	printf("===============================================================================\n");

	// n = 8
	const int dim2 = 8;
	double C2[dim2][dim2], Ct2[dim2][dim2], CtC2[dim2][dim2];
	double d_2[dim2], d1_2[dim2];
	int ipvt2[dim2];
	double cond2, delta2;
	double m2 = 1, n2 = 1;

	for (int i = 0; i < dim2; i++) {
		for (int j = 0; j < dim2; j++) {
			C2[i][j] = 1 / (m2 + n2 - 1);
			n2++;
		}
		n2 = 1; m2++;
	}

	m2 = 1; n2 = 1;
	for (int i = 0; i < dim2; i++) {
		d_2[i] = 0;
		for (int j = 0; j < dim2; j++) {
			d_2[i] = 1 / (m2 + n2 - 1) + d_2[i];
			n2++;
		}
		n2 = 1; m2++;
	}

	for (int i = 0; i < dim2; i++)
		for (int j = 0; j < dim2; j++)
			Ct2[i][j] = C2[j][i];

	m2 = 1; n2 = 1;
	for (int i = 0; i < dim2; i++) {
		d1_2[i] = 0;
		for (int j = 0; j < dim2; j++) {
			d1_2[i] = d1_2[i] + Ct2[i][j] * d_2[j];
			n2++;
		}
		n2 = 1; m2++;
	}

	for (int i = 0; i < dim2; i++) {
		for (int j = 0; j < dim2; j++) {
			CtC2[i][j] = 0;
			for (int k = 0; k < dim2; k++)
				CtC2[i][j] = CtC2[i][j] + Ct2[i][k] * C2[k][j];
		}
	}

	printf("      Matrix C                                      d\n");
	for (int i = 0; i < dim2; i++) {
		for (int j = 0; j < dim2; j++)
			printf("%1.2f  ", C2[i][j]);
		printf("|  %1.2f", d_2[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim2, *C2, &cond2, ipvt2);
	printf("cond = %e\n\n", cond2);
	Solve(dim2, *C2, d_2, ipvt2);
	int z2 = 1;
	for (int i = 0; i < dim2; i++) {
		printf("x1[%d] =  %1.2f", z2, d_2[i]);
		z2++;
		printf("\n");
	}
	printf("_______________________________________________________\n\n");

	printf("      Matrix CtC                                   Ct*d\n");
	for (int i = 0; i < dim2; i++) {
		for (int j = 0; j < dim2; j++)
			printf("%1.2f  ", CtC2[i][j]);
		printf("|  %1.2f", d1_2[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim2, *CtC2, &cond2, ipvt2);
	printf("cond = %e\n\n", cond2);
	Solve(dim2, *CtC2, d1_2, ipvt2);
	z2 = 1;
	for (int i = 0; i < dim2; i++) {
		printf("x2[%d] =  %1.2f", z2, d1_2[i]);
		z2++;
		printf("\n");
	}
	printf("_______________________________\n");

	double s1_2 = 0, s2_2 = 0;
	for (int i = 0; i < dim2; i++) {
		s1_2 = pow(d_2[i], 2) + s1_2;
		s2_2 = pow((d_2[i] - d1_2[i]), 2) + s2_2;
	}
	s1_2 = sqrt(s1_2);
	s2_2 = sqrt(s2_2);
	delta2 = s2_2 / s1_2;
	printf("delta = %e\n\n", delta2);
	printf("===============================================================================\n");

	// n = 10
	const int dim3 = 10;
	double C3[dim3][dim3], Ct3[dim3][dim3], CtC3[dim3][dim3];
	double d_3[dim3], d1_3[dim3];
	int ipvt3[dim3];
	double cond3, delta3;
	double m3 = 1, n3 = 1;

	for (int i = 0; i < dim3; i++) {
		for (int j = 0; j < dim3; j++) {
			C3[i][j] = 1 / (m3 + n3 - 1);
			n3++;
		}
		n3 = 1; m3++;
	}

	m3 = 1; n3 = 1;
	for (int i = 0; i < dim3; i++) {
		d_3[i] = 0;
		for (int j = 0; j < dim3; j++) {
			d_3[i] = 1 / (m3 + n3 - 1) + d_3[i];
			n3++;
		}
		n3 = 1; m3++;
	}

	for (int i = 0; i < dim3; i++)
		for (int j = 0; j < dim3; j++)
			Ct3[i][j] = C3[j][i];

	m3 = 1; n3 = 1;
	for (int i = 0; i < dim3; i++) {
		d1_3[i] = 0;
		for (int j = 0; j < dim3; j++) {
			d1_3[i] = d1_3[i] + Ct3[i][j] * d_3[j];
			n3++;
		}
		n3 = 1; m3++;
	}

	for (int i = 0; i < dim3; i++) {
		for (int j = 0; j < dim3; j++) {
			CtC3[i][j] = 0;
			for (int k = 0; k < dim3; k++)
				CtC3[i][j] = CtC3[i][j] + Ct3[i][k] * C3[k][j];
		}
	}

	printf("      Matrix C                                                  d\n");
	for (int i = 0; i < dim3; i++) {
		for (int j = 0; j < dim3; j++)
			printf("%1.2f  ", C3[i][j]);
		printf("|  %1.2f", d_3[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim3, *C3, &cond3, ipvt3);
	printf("cond = %e\n\n", cond3);
	Solve(dim3, *C3, d_3, ipvt3);
	int z3 = 1;
	for (int i = 0; i < dim3; i++) {
		printf("x1[%d] =  %1.2f", z3, d_3[i]);
		z3++;
		printf("\n");
	}
	printf("___________________________________________________________________\n\n");

	printf("      Matrix CtC                                               Ct*d\n");
	for (int i = 0; i < dim3; i++) {
		for (int j = 0; j < dim3; j++)
			printf("%1.2f  ", CtC3[i][j]);
		printf("|  %1.2f", d1_3[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim3, *CtC3, &cond3, ipvt3);
	printf("cond = %e\n\n", cond3);
	Solve(dim3, *CtC3, d1_3, ipvt3);
	z3 = 1;
	for (int i = 0; i < dim3; i++) {
		printf("x2[%d] =  %1.2f", z3, d1_3[i]);
		z3++;
		printf("\n");
	}
	printf("_______________________________\n");

	double s1_3 = 0, s2_3 = 0;
	for (int i = 0; i < dim3; i++) {
		s1_3 = pow(d_3[i], 2) + s1_3;
		s2_3 = pow((d_3[i] - d1_3[i]), 2) + s2_3;
	}
	s1_3 = sqrt(s1_3);
	s2_3 = sqrt(s2_3);
	delta3 = s2_3 / s1_3;
	printf("delta = %e\n\n", delta3);
	printf("===============================================================================\n");

	// n = 12
	const int dim4 = 12;
	double C4[dim4][dim4], Ct4[dim4][dim4], CtC4[dim4][dim4];
	double d_4[dim4], d1_4[dim4];
	int ipvt4[dim4];
	double cond4, delta4;
	double m4 = 1, n4 = 1;

	for (int i = 0; i < dim4; i++) {
		for (int j = 0; j < dim4; j++) {
			C4[i][j] = 1 / (m4 + n4 - 1);
			n4++;
		}
		n4 = 1; m4++;
	}

	m4 = 1; n4 = 1;
	for (int i = 0; i < dim4; i++) {
		d_4[i] = 0;
		for (int j = 0; j < dim4; j++) {
			d_4[i] = 1 / (m4 + n4 - 1) + d_4[i];
			n4++;
		}
		n4 = 1; m4++;
	}

	for (int i = 0; i < dim4; i++)
		for (int j = 0; j < dim4; j++)
			Ct4[i][j] = C4[j][i];

	m4 = 1; n4 = 1;
	for (int i = 0; i < dim4; i++) {
		d1_4[i] = 0;
		for (int j = 0; j < dim4; j++) {
			d1_4[i] = d1_4[i] + Ct4[i][j] * d_4[j];
			n4++;
		}
		n4 = 1; m4++;
	}

	for (int i = 0; i < dim4; i++) {
		for (int j = 0; j < dim4; j++) {
			CtC4[i][j] = 0;
			for (int k = 0; k < dim4; k++)
				CtC4[i][j] = CtC4[i][j] + Ct4[i][k] * C4[k][j];
		}
	}

	printf("      Matrix C                                                              d\n");
	for (int i = 0; i < dim4; i++) {
		for (int j = 0; j < dim4; j++)
			printf("%1.2f  ", C4[i][j]);
		printf("|  %1.2f", d_4[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim4, *C4, &cond4, ipvt4);
	printf("cond = %e\n\n", cond4);
	Solve(dim4, *C4, d_4, ipvt4);
	int z4 = 1;
	for (int i = 0; i < dim4; i++) {
		printf("x1[%d] =  %1.2f", z4, d_4[i]);
		z4++;
		printf("\n");
	}
	printf("_______________________________________________________________________________\n\n");

	printf("      Matrix CtC                                                           Ct*d\n");
	for (int i = 0; i < dim4; i++) {
		for (int j = 0; j < dim4; j++)
			printf("%1.2f  ", CtC4[i][j]);
		printf("|  %1.2f", d1_4[i]);
		printf("\n");
	}
	printf("\n");

	Decomp(dim4, *CtC4, &cond4, ipvt4);
	printf("cond = %e\n\n", cond4);
	Solve(dim4, *CtC4, d1_4, ipvt4);
	z4 = 1;
	for (int i = 0; i < dim4; i++) {
		printf("x2[%d] =  %1.2f", z4, d1_4[i]);
		z4++;
		printf("\n");
	}

	printf("_______________________________\n");

	double s1_4 = 0, s2_4 = 0;
	for (int i = 0; i < dim4; i++) {
		s1_4 = pow(d_4[i], 2) + s1_4;
		s2_4 = pow((d_4[i] - d1_4[i]), 2) + s2_4;
	}
	s1_4 = sqrt(s1_4);
	s2_4 = sqrt(s2_4);
	delta4 = s2_4 / s1_4;
	printf("delta = %e\n\n", delta4);

	system("pause");
}
