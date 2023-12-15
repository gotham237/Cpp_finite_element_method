#pragma once
#ifndef PRINT_DATA_TO_FILE_H 
#define PRINT_DATA_TO_FILE_H

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include "getDataFromFile.h"
#include "SoE.h"

using namespace std;

void createParaViewFile(Grid& grid, int i, SoE& soe) {
	std::string fileName = "file_" + std::to_string(i) + ".vtk";

	ofstream outputFile(fileName);
	if (outputFile.is_open()) {

		/*# vtk DataFile Version 2.0
			Unstructured Grid Example
			ASCII
			DATASET UNSTRUCTURED_GRID*/
		outputFile << "# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";
		outputFile << "POINTS " + std::to_string(grid.nN) + " float\n";

		for (int i = 0; i < grid.nN; i++) {
			outputFile << grid.tNode[i].x << " " << grid.tNode[i].y << " 0\n";
		}
		outputFile << "\n";

		outputFile << "CELLS " << grid.nE << " " << grid.nE * 5 << "\n";
		for (int i = 0; i < grid.nE; i++) {
			outputFile << "4 ";
			for (int j = 0; j < 4; j++) { // na razie 4 tyle co w ID jest
				outputFile << grid.tElem[i].ID[j] - 1 << " ";
			}
			outputFile << "\n";
		}
		outputFile << "\n";

		outputFile << "CELL_TYPES " << grid.nE << "\n";
		for (int i = 0; i < grid.nE; i++) {
			outputFile << 9 << "\n";
		}
		outputFile << "\n";

		outputFile << "POINT_DATA " << grid.nN << "\n";
		outputFile << "SCALARS Temp float 1" << "\n";
		outputFile << "LOOKUP_TABLE default " << "\n";

		for (int i = 0; i < grid.nN; i++) {
			outputFile << soe.newTemperatures[i] << "\n";
		}

		outputFile.close();
		std::cout << "Data has been written to the file." << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}
}

#endif
