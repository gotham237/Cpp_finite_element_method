#pragma once
#include "getDataFromFile.h";

double* Gauss(int wiersze, int kolumny, double** macierz);

struct SoE{
	double** H;
	double* P;
    double** C;
    double** gauss;
    double* modP;

    double* newTemperatures;

    void calculate(Grid grid, GlobalData& globalData) {
        const int nN = grid.nN;
        H = new double* [nN];
        P = new double[nN];
        C = new double* [nN];
        gauss = new double* [nN];
        newTemperatures = new double[nN];
        modP = new double[nN];

        for (int i = 0; i < nN; i++) {
            H[i] = new double[nN];
            C[i] = new double[nN];
            gauss[i] = new double[nN + 1];

            for (int j = 0; j < nN; j++) {
                H[i][j] = 0.0;
                C[i][j] = 0.0;
                gauss[i][j] = 0.0;
            }
            newTemperatures[i] = globalData.initialTemp;
            gauss[i][nN] = 0.0;
            P[i] = 0.0;
            modP[i] = 0.0;
        }

        const int nP = 4;
        for (int i = 0; i < grid.nE; i++) {
            for (int j = 0; j < nP; j++) {
                for (int k = 0; k < nP; k++) {
                    H[grid.tElem[i].ID[j] - 1][grid.tElem[i].ID[k] - 1] += grid.tElem[i].H[j][k];
                    C[grid.tElem[i].ID[j] - 1][grid.tElem[i].ID[k] - 1] += grid.tElem[i].C[j][k];
                }
                P[grid.tElem[i].ID[j] - 1] += grid.tElem[i].P[j];
                modP[grid.tElem[i].ID[j] - 1] += grid.tElem[i].P[j];
            }
        }
	}

    void calculateTemperature(Grid& grid) {
        for (int j = 0; j < grid.nN; j++) {
            for (int k = 0; k < grid.nN; k++) {
                gauss[j][k] = H[j][k];
            }
            gauss[j][grid.nN] = modP[j];
        }

        newTemperatures = Gauss(grid.nN, grid.nN + 1, gauss);
    }

    void calculateModP(Grid& grid, GlobalData& globalData) {
        double* tempP = new double[grid.nN];

        for (int i = 0; i < grid.nN; i++) {
            tempP[i] = 0.0;
            for (int j = 0; j < grid.nN; j++) {
                tempP[i] += (C[i][j] / globalData.simulationStepTime) * newTemperatures[j];
            }
        }

        for (int i = 0; i < grid.nN; i++) {
            modP[i] = P[i] + tempP[i];
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
        cout << "\n";

		cout << "P globalne: " << "\n";
		for (int i = 0; i < grid.nN; i++) {
			cout << modP[i] << "  ";
		}
		cout << "\n\n";

        cout << "C globalne: " << "\n";
        for (int i = 0; i < grid.nN; i++) {
            for (int j = 0; j < grid.nN; j++) {
                cout << C[i][j] << "  ";
            }
            cout << "\n";
        }
        cout << "\n";
	}
};

double* Gauss(int wiersze, int kolumny, double** macierz) {
    double mnoznik = 0.0;
    for (int k = 0; k < wiersze; k++) {
        for (int i = k; i < wiersze - 1; i++) {
            if (macierz[i][i] != 0.0) {
                mnoznik = macierz[i + 1][k] / macierz[k][k];
                for (int j = k; j < kolumny; j++) {
                    macierz[i + 1][j] -= macierz[k][j] * mnoznik;
                }
            }
            else {
                cout << "Nie mozna dzeilic przez 0" << endl;
                break;
            }
        }
    }

    double* tabx = new double[wiersze - 1];
    tabx[wiersze - 1] = macierz[wiersze - 1][wiersze] / macierz[wiersze - 1][wiersze - 1];
    double xn = 0.0;
    double sum = 0.0;

    for (int i = wiersze - 2; i >= 0; i--) {
        sum = 0.0;
        for (int k = i + 1; k < wiersze; k++) {
            sum += macierz[i][k] * tabx[k];
        }
        tabx[i] = (macierz[i][wiersze] - sum) / macierz[i][i];
    }
    cout << "Rozwiazania ukladu: " << endl;
    for (int i = 0; i < wiersze; i++) {
        cout << "x" << i << " = " << tabx[i] << endl;
    }

    return tabx;
}
