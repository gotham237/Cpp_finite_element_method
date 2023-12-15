#pragma once
#ifndef GET_DATA_FROM_FILE_H  // This is an include guard to prevent multiple inclusion
#define GET_DATA_FROM_FILE_H 

#define NUM_OF_POINTS 2

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include "elemUniw.h"
#include "getDataFromFile.h"


using namespace std;

struct GlobalData
{
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
};

struct Node
{
    double x, y;
    double temp;
    int BC;
};

struct Element
{
    int ID[4];
    double H[4][4];
    double Hbc[4][4];
    double P[4];
    double C[4][4];

    Element() {
        ID[0] = 0, ID[1] = 0, ID[2] = 0, ID[3] = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Hbc[i][j] = 0;
                C[i][j] = 0;
                H[i][j] = 0;
            }
            P[i] = 0;
        }
    }

    void createElem(Node* tNode, int integralPoints, GlobalData& globalData)
    {
        /*double tabX[4] = { 0.0 , 0.025, 0.025, 0.0 };
        double tabY[4] = { 0.0, 0.0, 0.025, 0.025 };*/
        double tabX[4] = { tNode[ID[0] - 1].x, tNode[ID[1] - 1].x, tNode[ID[2] - 1].x, tNode[ID[3] - 1].x };
        double tabY[4] = { tNode[ID[0] - 1].y, tNode[ID[1] - 1].y, tNode[ID[2] - 1].y, tNode[ID[3] - 1].y };
        /*double tabX[4] = { 0.0, 0.02, 0.09 , 0.0 };
        double tabY[4] = { 0.0, 0.0, 0.04, 0.02};*/
        const int nP = integralPoints;
        elemUniw elementUniw(nP);
        PW pw = metodaGaussa("2d", nP);

        double** dNdx = new double* [nP * nP];
        double** dNdy = new double* [nP * nP];

        for (int k = 0; k < nP * nP; k++) {
            dNdx[k] = new double[4];
            dNdy[k] = new double[4];
        }

        for (int m = 0; m < 4; m++) {

            for (int l = 0; l < 4; l++) {
                H[m][l] = 0;
            }
        }

        for (int i = 0; i < nP * nP; i++)
        {
            double dxdEta = elementUniw.arrEta[i][0] * tabX[0] + elementUniw.arrEta[i][1] * tabX[1] + elementUniw.arrEta[i][2] * tabX[2] + elementUniw.arrEta[i][3] * tabX[3];
            double dxdKsi = elementUniw.arrKsi[i][0] * tabX[0] + elementUniw.arrKsi[i][1] * tabX[1] + elementUniw.arrKsi[i][2] * tabX[2] + elementUniw.arrKsi[i][3] * tabX[3];
            double dydEta = elementUniw.arrEta[i][0] * tabY[0] + elementUniw.arrEta[i][1] * tabY[1] + elementUniw.arrEta[i][2] * tabY[2] + elementUniw.arrEta[i][3] * tabY[3];
            double dydKsi = elementUniw.arrKsi[i][0] * tabY[0] + elementUniw.arrKsi[i][1] * tabY[1] + elementUniw.arrKsi[i][2] * tabY[2] + elementUniw.arrKsi[i][3] * tabY[3];

            //double wyznacznik[2][2] = { {dxdKsi, dydKsi}
            //                         , {dxdEta, dydEta} };

            double wyznacznik[2][2] = { {dydEta, -dydKsi}
                                      , {-dxdEta, dxdKsi} };

            double det = dxdKsi * dydEta - (-dydKsi * (-dxdEta));
            //cout << "det: " << det << "\n";


            for (int j = 0; j < 4; j++)
            {
                dNdx[i][j] = 1.0 / det * (wyznacznik[0][0] * elementUniw.arrKsi[i][j] + wyznacznik[0][1] * elementUniw.arrEta[i][j]);
                dNdy[i][j] = 1.0 / det * (wyznacznik[1][0] * elementUniw.arrKsi[i][j] + wyznacznik[1][1] * elementUniw.arrEta[i][j]);
            }

            //H i C
            double *wages = new double[nP * nP];
            if (nP == 2) {
                wages[0] = 1;
                wages[1] = 1;
                wages[2] = 1;
                wages[3] = 1;
            }
            else if (nP == 3) {
                wages[0] = pw.W[0] * pw.W[0];
                wages[1] = pw.W[0] * pw.W[1];
                wages[2] = pw.W[0] * pw.W[2];

                wages[3] = pw.W[1] * pw.W[0];
                wages[4] = pw.W[1] * pw.W[1];
                wages[5] = pw.W[1] * pw.W[2];

                wages[6] = pw.W[2] * pw.W[0];
                wages[7] = pw.W[2] * pw.W[1];
                wages[8] = pw.W[2] * pw.W[2];
            }

            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    H[j][k] += wages[i] * globalData.conductivity * (dNdx[i][j] * dNdx[i][k] + dNdy[i][j] * dNdy[i][k]) * det;
                    C[j][k] += wages[i] * globalData.specificHeat * globalData.density * (elementUniw.tabN[i][j] * elementUniw.tabN[i][k]) * det;
                }
            }

        }

        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                Hbc[j][k] = 0.0;
            }
        }

        //HBc i wektor P
        for (int i = 0; i < 4; i++) {
            //sprawdzenie czy wystepuje warunek brzegowy
            if (tNode[ID[i] - 1].BC == 1 && tNode[ID[(i + 1) % 4] - 1].BC == 1) {
                double x1 = tNode[ID[i] - 1].x;
                double y1 = tNode[ID[i] - 1].y;
                double x2 = tNode[ID[(i + 1) % 4] - 1].x;
                double y2 = tNode[ID[(i + 1) % 4] - 1].y;

                double L = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
                double detJ = L / 2.0;

                //cout << "J: " << detJ << "\n";

                //Hbc
                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        for (int m = 0; m < nP; m++) {
                            Hbc[j][k] += globalData.alfa * pw.W[m] * (elementUniw.surface[i].N[m][j] * elementUniw.surface[i].N[m][k]) * detJ;
                        }
                    }
                }


                //P
                for (int m = 0; m < nP; m++) {
                    for (int j = 0; j < 4; j++) {
                        P[j] += globalData.alfa * pw.W[m] * elementUniw.surface[i].N[m][j] * globalData.tot * detJ;
                    }
                }

            }

            
        }

        // H lokalne z HBC i C/dT i P + C/dT * {t0}
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                H[j][k] += Hbc[j][k] + (C[j][k] / globalData.simulationStepTime);
                //P[j] += (C[j][k] / globalData.simulationStepTime) * globalData.initialTemp;
            }
        }

        for (int j = 0; j < nP * nP; j++) {
            delete[] dNdx[j];
            delete[] dNdy[j];
        }

        delete[] dNdx;
        delete[] dNdy;
    }

    void printElement(int index) {
        cout << "H" << index << ": " << endl;
        for (int m = 0; m < 4; m++) {
            for (int l = 0; l < 4; l++) {
                cout << H[m][l] << "  ";
            }
            cout << "\n";
        }
        cout << "\n\n";

        cout << "Hbc" << index << ": \n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << Hbc[i][j] << "  ";
            }
            cout << "\n";
        }
        cout << "\n\n";

        cout << "P" << index << ": \n";
        for (int i = 0; i < 4; i++) {
            cout << P[i] << "  ";
        }
        cout << "\n\n";

        cout << "C" << index << ": \n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << C[i][j] << "  ";
            }
            cout << "\n";
        }
        cout << "\n\n";
    }
};

struct Grid
{
    int nN;
    int nE;
    Node* tNode;
    Element* tElem;
};

void getData(const string& filename, GlobalData& globalData, Grid& grid) 
{
    ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    string key;
    int value;

    while (inputFile >> key) {
        if (key == "SimulationTime")
            inputFile >> globalData.simulationTime;
        else if (key == "SimulationStepTime")
            inputFile >> globalData.simulationStepTime;
        else if (key == "Conductivity")
            inputFile >> globalData.conductivity;
        else if (key == "Alfa")
            inputFile >> globalData.alfa;
        else if (key == "Tot")
            inputFile >> globalData.tot;
        else if (key == "InitialTemp")
            inputFile >> globalData.initialTemp;
        else if (key == "Density")
            inputFile >> globalData.density;
        else if (key == "SpecificHeat")
            inputFile >> globalData.specificHeat;
        else if (key == "Nodes") {
            inputFile >> key >> grid.nN;
            grid.tNode = new Node[grid.nN];
        }
        else if (key == "Elements")
        {
            inputFile >> key >> grid.nE;
            grid.tElem = new Element[grid.nE];
        }
        else if (key == "*Node")
        {
            for (int i = 0; i < grid.nN; i++)
            {   
                string t;
                inputFile >> key >> t >> grid.tNode[i].y;
                t.pop_back();
                grid.tNode[i].x = stod(t);
                grid.tNode[i].temp = globalData.initialTemp;
            }
        }
        else if (key == "*BC")
        {
            string t;
            int i;
            while (inputFile >> t) {
                if (t[t.length() - 1] == ',') {
                    t.pop_back();
                }
                i = stoi(t);
                grid.tNode[i - 1].BC = 1;
            }
            for (int i = 0; i < grid.nN; i++) {
                if (grid.tNode[i].BC != 1) {
                    grid.tNode[i].BC = 0;
                }
            }
        }
    }
    inputFile.close();

    ifstream inputFile2(filename);

    if (!inputFile2.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    while (inputFile2 >> key) {
        if (key == "*Element,")
        {
            string one, two, three, four;
            inputFile2 >> key;
            for (int i = 0; i < grid.nE; i++)
            {
                inputFile2 >> key >> one >> two >> three >> four;
                grid.tElem[i].ID[0] = stod(one);
                grid.tElem[i].ID[1] = stod(two);
                grid.tElem[i].ID[2] = stod(three);
                grid.tElem[i].ID[3] = stod(four);

                grid.tElem[i].createElem(grid.tNode, NUM_OF_POINTS, globalData);
                //grid.tElem[i].printElement();
            }
        }
    }

    inputFile2.close();
}

void showGlobalData(GlobalData& globalData)
{
    cout << "SimulationTime: " << globalData.simulationTime << endl;
    cout << "SimulationStepTime: " << globalData.simulationStepTime << endl;
    cout << "Conductivity: " << globalData.conductivity << endl;
    cout << "Alfa: " << globalData.alfa << endl;
    cout << "Tot: " << globalData.tot << endl;
    cout << "InitialTemp: " << globalData.initialTemp << endl;
    cout << "Density: " << globalData.density << endl;
    cout << "SpecificHeat: " << globalData.specificHeat << endl << endl;
}

void showGridData(Grid& grid)
{
    cout << "Number of nodes: " << grid.nN << "\n";
    cout << "Number of elements: " << grid.nE << "\n";
    cout << "Nodes: " << "\n";
    for (int i = 0; i < grid.nN; i++)
    {
        cout << grid.tNode[i].x << "     " << grid.tNode[i].y << "    " << "\n";
    }
    cout << "\nBC: " << "\n";
    for (int i = 0; i < grid.nN; i++) {
        cout << grid.tNode[i].BC << "  ";
    }
    cout << "\n\n";
    cout << "Elements: " << "\n";
    for (int i = 0; i < grid.nE; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << grid.tElem[i].ID[j] << "  ";
        }
        cout << "\n";
    }
    cout << "\n";

    for (int i = 0; i < grid.nE; i++) {
        grid.tElem[i].printElement(i + 1);
    }
}



#endif