#pragma once
#ifndef ELEM_UNIW_H
#define ELEM_UNIW_H

#include <iostream>
#include <string>
#include "metodaGaussa.h"

struct Surface {
    double** N;
};

struct elemUniw {
    int nP;
    //pochodne po ksi i eta
    double** arrKsi;
    double** arrEta;
    Surface surface[4];

    elemUniw(int integralPoints)
    {
        nP = integralPoints;
        double N = nP * nP;

        arrKsi = new double* [N];
        arrEta = new double* [N];

        for (int i = 0; i < N; i++) {
            arrKsi[i] = new double[4];
            arrEta[i] = new double[4];
        }
        //points and wages
        PW pw = metodaGaussa("2d", nP);
        //double* points = metodaGaussa("2d", nP);

        for (int i = 0; i < N; i++)
        {

            int* iE = new int[N];
            int* iKsi = new int[N];

            if (nP == 2) {
                iE[0] = 0; iE[1] = 0; iE[2] = 1; iE[3] = 1;
                iKsi[0] = 0; iKsi[1] = 1; iKsi[2] = 0; iKsi[3] = 1;
            }
            else if (nP == 3) {
                iE[0] = 0; iE[1] = 1; iE[2] = 2;
                iE[3] = 0; iE[4] = 1; iE[5] = 2;
                iE[6] = 0; iE[7] = 1; iE[8] = 2;

                iKsi[0] = 0; iKsi[1] = 0; iKsi[2] = 0;
                iKsi[3] = 1; iKsi[4] = 1; iKsi[5] = 1;
                iKsi[6] = 2; iKsi[7] = 2; iKsi[8] = 2;
            }
            else if (nP == 4) {
                iE[0] = 0; iE[1] = 0; iE[2] = 0; iE[3] = 0;
                iE[4] = 1; iE[5] = 1; iE[6] = 1; iE[7] = 1;
                iE[8] = 2; iE[9] = 2; iE[10] = 2; iE[11] = 2;
                iE[12] = 3; iE[13] = 3; iE[14] = 3; iE[15] = 3;

                iKsi[0] = 0; iKsi[1] = 1; iKsi[2] = 2; iKsi[3] = 3;
                iKsi[4] = 0; iKsi[5] = 1; iKsi[6] = 2; iKsi[7] = 3;
                iKsi[8] = 0; iKsi[9] = 1; iKsi[10] = 2; iKsi[11] = 3;
                iKsi[12] = 0; iKsi[13] = 1; iKsi[14] = 2; iKsi[15] = 3;
            }

            arrKsi[i][0] = -0.25 * (1 - pw.P[iE[i]]);
            arrKsi[i][1] = 0.25 * (1 - pw.P[iE[i]]);
            arrKsi[i][2] = 0.25 * (1 + pw.P[iE[i]]);
            arrKsi[i][3] = -0.25 * (1 + pw.P[iE[i]]);

            arrEta[i][0] = -0.25 * (1 - pw.P[iKsi[i]]);
            arrEta[i][1] = -0.25 * (1 + pw.P[iKsi[i]]);
            arrEta[i][2] = 0.25 * (1 + pw.P[iKsi[i]]);
            arrEta[i][3] = 0.25 * (1 - pw.P[iKsi[i]]);

            

            delete[] iE;
            delete[] iKsi;

        }

        double Nodes[4][2][2] = {
                    {{-1.0 / sqrt(3.0), -1.0}, {1.0 / sqrt(3.0), -1.0}},
                    {{1.0, -1.0 / sqrt(3.0)},  {1.0, 1.0 / sqrt(3.0)}},
                    {{1.0 / sqrt(3.0), 1.0},   {-1.0 / sqrt(3.0), 1.0}},
                    {{-1.0 , 1.0 / sqrt(3.0)},  {-1.0, -1.0 / sqrt(3.0)}}
        };

        double Nodes3[4][3][2] = {
                    {{-sqrt(3.0 / 5.0), -1.0}, {0.0, -1.0}, {sqrt(3.0 / 5.0), -1.0}},
                    {{1.0, -sqrt(3.0 / 5.0)}, {1.0, 0.0}, {1.0, sqrt(3.0 / 5.0)}},
                    {{sqrt(3.0 / 5.0), 1}, {0.0, 1.0}, {-sqrt(3.0 / 5.0), 1}},
                    {{-1.0, sqrt(3.0 / 5.0)}, {-1.0, 0.0}, {-1.0, -sqrt(3.0 / 5.0)}}
        };

        //double edges[4] = { -1, 1, 1, -1 };
        for (int i = 0; i < 4; i++) {
            surface[i].N = new double* [nP];

            for (int j = 0; j < nP; j++) {
                surface[i].N[j] = new double[4];

                if (nP == 2) {
                    surface[i].N[j][0] = 0.25 * (1 - Nodes[i][j][0]) * (1 - Nodes[i][j][1]);
                    surface[i].N[j][1] = 0.25 * (1 + Nodes[i][j][0]) * (1 - Nodes[i][j][1]);
                    surface[i].N[j][2] = 0.25 * (1 + Nodes[i][j][0]) * (1 + Nodes[i][j][1]);
                    surface[i].N[j][3] = 0.25 * (1 - Nodes[i][j][0]) * (1 + Nodes[i][j][1]);
                }
                else if (nP == 3) {
                    surface[i].N[j][0] = 0.25 * (1 - Nodes3[i][j][0]) * (1 - Nodes3[i][j][1]);
                    surface[i].N[j][1] = 0.25 * (1 + Nodes3[i][j][0]) * (1 - Nodes3[i][j][1]);
                    surface[i].N[j][2] = 0.25 * (1 + Nodes3[i][j][0]) * (1 + Nodes3[i][j][1]);
                    surface[i].N[j][3] = 0.25 * (1 - Nodes3[i][j][0]) * (1 + Nodes3[i][j][1]);
                }
            }
        }

    }

    void printSurface() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 4; k++) {
                    cout << surface[i].N[j][k] << "   ";
                }
                cout << endl;
            }
            cout << "\n\n";
        }
    }


    void printElemUniw()
    {
        std::cout << "KSI:\n";
        for (int i = 0; i < nP * nP; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << arrKsi[i][j] << " ";
            }
            std::cout << "\n";
        }

        std::cout << "ETA:\n";
        for (int i = 0; i < nP * nP; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << arrEta[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

    ~elemUniw()
    {
        for (int i = 0; i < nP * nP; i++)
        {
            delete[] arrKsi[i];
            delete[] arrEta[i];
        }

        delete[] arrKsi;
        delete[] arrEta;
    }
};
#endif