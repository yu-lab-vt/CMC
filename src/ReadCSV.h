#pragma once


#ifndef READCSV_H
#define READCSV_H



#include "Global.h"




uint* ReadCSV_UINT(char* fileName, uint nHeaderlines); // UINT (1D)
double* ReadCSV_DOUBLE(char* fileName, uint nHeaderlines, uint* p_nRow); // double (1D)
uint** ReadCSV_UINT_2D(char* fileName, uint nHeaderlines); // UINT (2D)
void readBin(double** buff, char* fileName, uint* size);
void readBin(uint** buff, char* fileName, uint* size);
void readBin(bool** buff, char* fileName, unsigned long long int* size);
#endif



