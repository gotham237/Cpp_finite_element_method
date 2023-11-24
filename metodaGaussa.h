#pragma once
#ifndef METODA_GAUSSA_H
#define METODA_GAUSSA_H 
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

struct PW
{
	double* P;
	double* W;

	PW(int nP)
	{
		P = new double[nP];
		W = new double[nP];
	}
};

double f1(double x)
{
	return 5 * x * x + 3 * x + 6;
}

double f2(double x, double y)
{
	return 5 * x * x * y * y + 3 * x * y + 6;
}

PW metodaGaussa(string wymiar, int nP)
{
	PW pw = PW(nP);
	if (nP == 2)
	{
		pw.P[0] = -sqrt(1.0 / 3.0);
		pw.P[1] = sqrt(1.0 / 3.0);
		pw.W[0] = 1.0;
		pw.W[1] = 1.0;
	}
	else if (nP == 3)
	{
		pw.P[0] = sqrt(3.0 / 5.0);
		pw.P[1] = 0.0;
		pw.P[2] = -1.0 * sqrt(3.0 / 5.0);
		pw.W[0] = 5.0 / 9.0;
		pw.W[1] = 8.0 / 9.0;
		pw.W[2] = 5.0 / 9.0;
	}
	else if (nP == 4)
	{
		pw.P[0] = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.P[1] = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.P[2] = -1 * sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.P[3] = -1 * sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
		pw.W[0] = (18 + sqrt(30.0)) / 36.0;
		pw.W[3] = (18 - sqrt(30.0)) / 36.0;
		pw.W[1] = (18 - sqrt(30.0)) / 36.0;
		pw.W[2] = (18 + sqrt(30.0)) / 36.0;
	}

	double suma = 0.0;
	if (wymiar == "1d")
	{
		for (int i = 0; i < nP; i++)
		{
			suma += pw.W[i] * f1(pw.P[i]);
		}
	}
	else if (wymiar == "2d")
	{
		for (int i = 0; i < nP; i++)
		{
			for (int j = 0; j < nP; j++)
			{
				suma += pw.W[i] * pw.W[j] * f2(pw.P[i], pw.P[j]);
			}
		}
	}
	return pw;
}

#endif
