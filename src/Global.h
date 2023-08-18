#pragma once



#ifndef Global_H
#define Global_H

#include <iostream>
#include<fstream>
#include<string>
using namespace std;

#include "Defines.h"




uint** CALLOC_2D_UINT(uint nRow, uint nCol);






/************* Class *******************/

/************* Class *******************/

class Index_HD {  // Index of high dimensions
//private:

public:
    uint nDim;
    uint M[MAX_DIM]; //length of each dimension
    uint idx[MAX_DIM] = { 0 };
    uint iExcess = 0;

    // Constructor function
    Index_HD(const uint i_nDim, const uint* i_M)
    {
        nDim = i_nDim;
        for (uint iDim = 0; iDim < nDim; iDim++) {
            if (i_M[iDim] == 0) {
                std::cout << "Dimension with length equal to 0 is not allowed!\n";
                exit(0);
            }

            M[iDim] = i_M[iDim];
        }
    }

    // Operater ++
    uint* operator ++ (int) {
        if (iExcess == 1)
            return idx;
        iExcess = 1;
        for (int iDim = 0; iDim < nDim; iDim++) {
            if (++idx[iDim] < M[iDim]) {
                iExcess = 0;
                break;
            }
            else {
                idx[iDim] = 0;
            }

        }
        return idx;
    }


    void reset() {
        for (uint iDim = 0; iDim < nDim; iDim++)
            idx[iDim] = 0;
    }

    uint getCurState() { return iExcess; }

    uint* getCurIdx() { return idx; }

};



class R_DIM {

public:
    double* r[MAX_DIM] = { NULL };
    double* r_pre[MAX_DIM] = { NULL };
    double* delta_r[MAX_DIM] = { NULL };
    uint nDim = 0;
    uint M[MAX_DIM] = { 0 }; //length of each dimension


    // Constructor function
    R_DIM(const uint i_nDim, const uint* i_M) {
        uint iDim = 0;
        nDim = i_nDim;
        for (iDim = 0; iDim < nDim; iDim++) {
            M[iDim] = i_M[iDim];
            r[iDim] = CALLOC(double, M[iDim]);
            r_pre[iDim] = CALLOC(double, M[iDim]);
            delta_r[iDim] = CALLOC(double, M[iDim]);
            for (uint j = 0; j < M[iDim]; j++) { // Initialization
                r[iDim][j] = 1.0;
                r_pre[iDim][j] = 1.0;
                delta_r[iDim][j] = 0.0;
            }
        }
    };

    R_DIM() {};


    // other functions
    void update(uint iDim, uint idx, double val) {
        r_pre[iDim][idx] = r[iDim][idx];
        delta_r[iDim][idx] = val - r[iDim][idx];
        r[iDim][idx] = val;
    }


    R_DIM& reorderElement(uint** idx) {
        // Reorder
        uint iDim = 0, j = 0;
        double* r_tmp = NULL;
        double* r_pre_tmp = NULL;
        double* delta_r_tmp = NULL;
        for (iDim = 0; iDim < this->nDim; iDim++) {
            r_tmp = CALLOC(double, this->M[iDim]);
            r_pre_tmp = CALLOC(double, this->M[iDim]);
            delta_r_tmp = CALLOC(double, this->M[iDim]);
            for (j = 0; j < M[iDim]; j++) { // tmp save
                r_tmp[j] = this->r[iDim][j];
                r_pre_tmp[j] = this->r_pre[iDim][j];
                delta_r_tmp[j] = this->delta_r[iDim][j];
            }
            for (j = 0; j < M[iDim]; j++) { // reorder
                this->r[iDim][j] = r_tmp[idx[iDim][j]];
                this->r_pre[iDim][j] = r_pre_tmp[idx[iDim][j]];
                this->delta_r[iDim][j] = delta_r_tmp[idx[iDim][j]];
            }
            free(r_tmp);
            free(r_pre_tmp);
            free(delta_r_tmp);
        }
        return *this;
    }

    void reset() {
        uint iDim = 0;
        for (iDim = 0; iDim < nDim; iDim++) {
            for (uint j = 0; j < M[iDim]; j++) { // Initialization
                r[iDim][j] = 1.0;
                r_pre[iDim][j] = 1.0;
                delta_r[iDim][j] = 0.0;
            }
        }
    }




    R_DIM& operator= (const R_DIM& r2) {
        uint iDim = 0;

        // FREE previous memory
        for (iDim = 0; iDim < nDim; iDim++) {
            FREE(r[iDim]);
            FREE(r_pre[iDim]);
            FREE(delta_r[iDim]);
        }
            

        // Copy
        nDim = r2.nDim;
        for (iDim = 0; iDim < nDim; iDim++) {
            M[iDim] = r2.M[iDim];
            r[iDim] = CALLOC(double, M[iDim]);
            r_pre[iDim] = CALLOC(double, M[iDim]);
            delta_r[iDim] = CALLOC(double, M[iDim]);
            for (uint j = 0; j < M[iDim]; j++) { // Initialization
                r[iDim][j] = r2.r[iDim][j];
                r_pre[iDim][j] = r2.r_pre[iDim][j];
                delta_r[iDim][j] = r2.delta_r[iDim][j];
            }
        }
        return *this;
    }


    double operator- (const R_DIM& r2) {
        if (nDim != r2.nDim) {
            std::cout << "ERROR! The numbers two Rs' dimensions are not the same.\n";
            exit(0);
        }

        double diff = 0.0;
        double maxDiff = 0.0;
        for (uint iDim = 0; iDim < nDim; iDim++) {
            if (M[iDim] != r2.M[iDim]) {
                std::cout << "ERROR! The lengths of two Rs' dimension are not the same.\n";
                exit(0);
            }

            for (uint j = 0; j < M[iDim]; j++) {
                diff = abs(r[iDim][j] - r2.r[iDim][j]);
                if (diff > maxDiff)
                    maxDiff = diff;
            }
        }
        return maxDiff;
    }


    // Destructor function
    ~R_DIM() {
        for (uint iDim = 0; iDim < nDim; iDim++) {
            FREE(r[iDim]);
            FREE(r_pre[iDim]);
            FREE(delta_r[iDim]);
        }
    }

};




#endif




