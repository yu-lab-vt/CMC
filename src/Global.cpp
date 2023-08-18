

#include "Global.h"


uint** CALLOC_2D_UINT(uint nRow, uint nCol) {

	uint len = sizeof(uint*) * nRow + sizeof(uint) * nCol * nRow;
	uint** matrix_2D = (uint**)malloc(len);

	// ptr is now pointing to the first element in of 2D array
	uint* ptr = (uint*)(matrix_2D + nRow);

	// for loop to point rows pointer to appropriate location in 2D array
	for (uint i = 0; i < nRow; i++)
		matrix_2D[i] = (ptr + nCol * i);

	return matrix_2D;

}







