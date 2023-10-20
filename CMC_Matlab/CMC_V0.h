#pragma once


#ifndef CMC_V0_H
#define CMC_V0_H


#include "Global.h"
#include "ReadCSV.h"

void CMC(R_DIM* r_out, uint nDim, double** Margin, uint* M, uint len_Y, unsigned long long int len_Miss_);

#endif

