#include "getDataFromFile.h"
#include "elemUniw.h"
#include "metodaGaussa.h"
#include "SoE.h"

int main()
{
	cout << setprecision(10);

	GlobalData globalData;
	Grid grid = Grid();

	getData("Test2_4_4_MixGrid.txt", globalData, grid);
	showGlobalData(globalData);
	showGridData(grid);
	SoE soe = SoE();
	soe.calculate(grid);
	soe.printSoE(grid);
	//grid.tElem->printElement();
	
	/*elemUniw elementUniwersalny(2);
	elementUniwersalny.printSurface();*/
	

	return 0;
}