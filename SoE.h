#pragma once
#include "getDataFromFile.h";

struct SoE{
	double** H;
	double* P;

	void calculate(Grid grid) {
		const int nN = grid.nN;
		H = new double* [nN];
		P = new double[nN];
		for (int i = 0; i < nN; i++) {
			H[i] = new double[nN];
			for (int j = 0; j < nN; j++) {
				H[i][j] = 0;
			}
			P[i] = 0;
		}

		const int nP = 4;
		for (int i = 0; i < grid.nE; i++) {
			for (int j = 0; j < nP; j++) {
				for (int k = 0; k < nP; k++) {
					H[grid.tElem[i].ID[j] - 1][grid.tElem[i].ID[k] - 1] += grid.tElem[i].H[j][k];
				}
				P[grid.tElem[i].ID[j] - 1] += grid.tElem[i].P[j];
			}
		}
	}

	void printSoE(Grid grid) {
		cout << "H globalne: " << "\n";
		for (int i = 0; i < grid.nN; i++) {
			for (int j = 0; j < grid.nN; j++) {
				cout << H[i][j] << "  ";
			}
			cout << "\n";
		}

		cout << "P globalne: " << "\n";
		for (int i = 0; i < grid.nN; i++) {
			cout << P[i] << "  ";
		}
		cout << "\n";
	}
};
