#include "getDataFromFile.h"
#include "elemUniw.h"
#include "metodaGaussa.h"
#include "SoE.h"
#include "printDataToFile.h"

int main()
{
	cout << setprecision(10);

	GlobalData globalData;
	Grid grid = Grid();
	string file1 = "Test1_4_4.txt";
	string file2 = "Test2_4_4_MixGrid.txt";
	string file3 = "Test3_31_31_kwadrat.txt";

	getData(file2, globalData, grid);
	showGlobalData(globalData);
	showGridData(grid);

	SoE soe = SoE();
	soe.calculate(grid, globalData);

	int iterations = globalData.simulationTime / globalData.simulationStepTime;

	createParaViewFile(grid, 0, soe);
	for (int i = 1; i <= iterations; i++) {
		soe.calculateModP(grid, globalData);
		soe.calculateTemperature(grid);

		createParaViewFile(grid, i, soe);
	}

	//soe.printSoE(grid);

	return 0;
}