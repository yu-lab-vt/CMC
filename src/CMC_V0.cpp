// CMC_V0.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#define _FILE_OFFSET_BITS 64

#include <iostream>
//#include <math.h>       /* exp */
#include "CMC_V0.h"
#include "r_solver.cpp"
#include <chrono>// for running time
#include <ctime>  
#include <algorithm> // for sort
//#include <vector>  // for vector

#include <Rcpp.h>
using namespace Rcpp;


void solveEquationSet(R_DIM* p_r, double** Margin, const uint nDim, const uint* M);
void solveEquationSet_1Dim(int iIter_outer, R_DIM* p_r, const uint curDim, double** Margin, const uint nDim, const uint* M);
double solve_r(int Direction, bool iSameDim, double r_pre, double** r_Part, const uint nDim_Part, const uint* M_Part, const double m, const uint M_cur, double errT_r = 3.75e-10);
void f_df(double* f, double* df, double r_try, double** r_Part, const uint nDim_Part, const uint* M_Part, const double m);


uint iDebug = 0;
uint* Y = NULL;
uint flag_Y_unique = 0;
unsigned long long int idx_Y_based[MAX_DIM + 1] = { 0 };
bool* Miss_ = NULL;
uint flag_Miss = 0;



// [[Rcpp::export]]
List CMC(int nDim_input, NumericVector M_input, NumericVector Margin_input, NumericVector Y_input=1, NumericVector item_used_input=0)
{


    /* ****************************************************************************************
    // Get Input
    ***************************************************************************************** */
    //uint nDim = 3;
    uint nDim = nDim_input;
    uint iDim = 0, i = 0, j = 0, k = 0;
    uint nElements = 1, nVar = 0;



    /* Read tensor information (dimension, number of entries, etc) */
    uint M[MAX_DIM] = { 0 }; //Length of each dimension
    for (iDim = 0; iDim < nDim; iDim++)
        M[iDim] = M_input[iDim];
    for (nElements = 1, iDim = 0; iDim < nDim; iDim++) //Number of elements(entries)
        nElements = nElements * M[iDim];
    for (nVar = 0, iDim = 0; iDim < nDim; iDim++) //Number of variables (total length)
        nVar += M[iDim];


    /* Initilaze the output */
    List Result;
    NumericVector R_out(nVar);
    for (k = 0, iDim = 0; iDim < nDim; iDim++) {
        for (i = 0; i < M[iDim]; i++) {
            R_out[k] = -1;
            k++;
        }
    }
    Result["rwu"] = R_out;
    Result["ErrorFlag"] = 1;
    Result["ErrorText"] = "No error if ErrorFlag = 0";
    //   R_out.~NumericVector();
    //   return(Result);


    /* Read marginal information */
    double* Margin[MAX_DIM];  // Marginal total counts
    for (k = 0, iDim = 0; iDim < nDim; iDim++) {
        Margin[iDim] = CALLOC(double, M[iDim]);
        for (i = 0; i < M[iDim]; i++) {
            Margin[iDim][i] = Margin_input[k];
            k++;
        }
    }


    /* Read marginal information */
    /*
    char fileName_m_base[] = "L:/CMC_C_R_code/data/m0.bin";
    double* Margin[MAX_DIM];  // Marginal total counts
    uint M[MAX_DIM] = { 0 }; //Length of each dimension
    size_t len = strlen(fileName_m_base);
    for (iDim = 0; iDim < nDim; iDim++) {
        readBin(&Margin[iDim], fileName_m_base, &M[iDim]);
        fileName_m_base[len - 5] ++;
    }
    */




    /* Read Y */
    uint Y0 = 1;
    unsigned long long int  len_Y = Y_input.size();
    //char fileName_Y[] = "L:/CMC_C_R_code/data/Y_tmp.bin";
    //readBin(&Y, fileName_Y, &len_Y);
    Y = CALLOC(uint, len_Y);
    for (unsigned long long int i_tmp = 0; i_tmp < len_Y; i_tmp++)
        Y[i_tmp] = Y_input[i_tmp];

    // Decide which type of Y it is
    if (len_Y == 1) {
        Y0 = Y[0];
        if (Y[0] == 1)
            flag_Y_unique = 2; // Y_ij: all 1
        else
            flag_Y_unique = 1; // Y_ij: all the same, but not equal to 1
    }
    else {
        flag_Y_unique = 0; // Y_ij with different values for all ij
    }
    //flag_Y_unique = 2;


    /* Read Missing inf */
    unsigned long long int len_Miss_ = 0;
    /*
    flag_Miss = 0;
    if (flag_Miss) {
        char fileName_Miss_[] = "L:/CMC_C_R_code/data/Miss__pro.bin";
        readBin(&Miss_, fileName_Miss_, &len_Miss_);
    }*/
    len_Miss_ = item_used_input.size();
    Miss_ = CALLOC(bool, len_Miss_);
    for (unsigned long long int i_tmp = 0; i_tmp < len_Miss_; i_tmp++)
        Miss_[i_tmp] = item_used_input[i_tmp];
    if (len_Miss_ == 1)
        flag_Miss = 0;
    else
        flag_Miss = 1;


    /* Display the dimension information */
    Rcpp::Rcout << "Number of dimensions:" << nDim << std::endl;
    for (iDim = 0; iDim < nDim-1; iDim++)
        Rcpp::Rcout << M[iDim] << " X ";
    Rcpp::Rcout << M[iDim] << std::endl;
    Rcpp::Rcout << "The type of Y is: type " << flag_Y_unique << std::endl;
    if (flag_Miss==0)
        Rcpp::Rcout << "The data has no missing values" << std::endl;
    else
        Rcpp::Rcout << "The data contain some missing values" << std::endl;

    //return Result;

    /* Check & Preprocessing */
    auto start = std::chrono::system_clock::now(); // Time start

    // Check. The sum of margins in each dimension should be equal
    double sumMargin[MAX_DIM] = { 0 };
    for (iDim = 0; iDim < nDim; iDim++) {
        for (j = 0; j < M[iDim]; j++) {
            sumMargin[iDim] += Margin[iDim][j];
        }
    }
    for (i = 0; i < nDim; i++) {
        if (abs((sumMargin[i] - sumMargin[0]) / sumMargin[0]) > 0.0001) {
            Rcpp::Rcout << "ERROR! The sum of margins in each dimension should be equal." << std::endl;
            Result["ErrorText"] = "ERROR! The sum of margins in each dimension should be equal.";
            R_out.~NumericVector();
            return Result;
        }
    }


    // Check the length of Y & Miss_
    unsigned long long int nEntries = 1;
    for (iDim = 0; iDim < nDim; iDim++)
        nEntries = nEntries * M[iDim];
    if (flag_Y_unique == 0 && nEntries != len_Y) {
        Rcpp::Rcout << "ERROR! The length of Y is not consistent with the inputted margin." << std::endl;
        Result["ErrorText"] = "ERROR! The length of Y is not consistent with the inputted margin.";
        R_out.~NumericVector();
        return Result;
    }
    if (flag_Miss && nEntries != len_Miss_) {
        Rcpp::Rcout << "ERROR! The length of Miss_ is not consistent with the inputted margin." << std::endl;
        Result["ErrorText"] = "ERROR! The length of Miss_ is not consistent with the inputted margin.";
        R_out.~NumericVector();
        return Result;
    }



    /* Sorting marginal totals(for speed up in latter step) */
    double* Margin_s[MAX_DIM]; //Sort the marginal totals
    uint* Index_invS[MAX_DIM]; //Margin[iDim]=Margin_s[iDim][Index_invS[iDim]]
    uint* Index_s[MAX_DIM]; //Margin_s[iDim]=Margin[iDim][Index_s[iDim]]
    for (iDim = 0; iDim < nDim; iDim++) {
        Margin_s[iDim] = CALLOC(double, M[iDim]);
        Index_invS[iDim] = CALLOC(uint, M[iDim]);
        Index_s[iDim] = CALLOC(uint, M[iDim]);
        for (uint i = 0; i < M[iDim]; i++) {
            Margin_s[iDim][i] = Margin[iDim][i];
            Index_s[iDim][i] = i;
            Index_invS[iDim][i] = i;
        }
        //*****************
        sort(Index_s[iDim], Index_s[iDim] + M[iDim], \
            [&](int ii, int jj) {return Margin_s[iDim][ii] < Margin_s[iDim][jj]; });
        sort(Index_invS[iDim], Index_invS[iDim] + M[iDim], \
            [&](int ii, int jj) {return Index_s[iDim][ii] < Index_s[iDim][jj]; });
        sort(Margin_s[iDim], Margin_s[iDim] + M[iDim]);
    }

    /* Sort Y according to Index_s*/
    if (flag_Y_unique == 0 || flag_Miss) {
        unsigned long long int* idx_base = CALLOC(unsigned long long int, nDim + 1);
        idx_base[0] = 1;
        for (iDim = 0; iDim < nDim - 1; iDim++) {
            idx_base[iDim + 1] = idx_base[iDim] * M[iDim];
        }
        Index_HD idx_Y(nDim, M);
        idx_Y.reset();
        unsigned long long int i_tmp = 0, idx_tmp = 0;
        uint* Y_s = CALLOC(uint, len_Y);
        bool* Miss__s = CALLOC(bool, len_Miss_);
        for (i_tmp = 0; i_tmp < nEntries; i_tmp++) {
            idx_tmp = 0;
            for (iDim = 0; iDim < nDim; iDim++) {
                idx_tmp += idx_base[iDim] * Index_s[iDim][idx_Y.idx[iDim]];
            }
            if (flag_Y_unique == 0)
                Y_s[i_tmp] = Y[idx_tmp];
            if (flag_Miss)
                Miss__s[i_tmp] = Miss_[idx_tmp];
            idx_Y++;
        }
        if (flag_Y_unique == 0)
            for (i_tmp = 0; i_tmp < len_Y; i_tmp++)
                Y[i_tmp] = Y_s[i_tmp];
        if (flag_Miss)
            for (i_tmp = 0; i_tmp < len_Miss_; i_tmp++)
                Miss_[i_tmp] = Miss__s[i_tmp];

        FREE(Y_s);
        FREE(Miss__s);
        FREE(idx_base);
    }


    /* Divide the margins by Y_ij, if Y_ij are all the same but not equal to 1 */
    if (flag_Y_unique == 1) {
        for (uint iDim = 0; iDim < nDim; iDim++) {
            for (j = 0; j < M[iDim]; j++) {
                Margin_s[iDim][j] = Margin_s[iDim][j] / Y0;
            }
        }
    }



    /* ****************************************************************************************
    // CMC
    ***************************************************************************************** */
    R_DIM r0(nDim, M);
    //solveEquationSet;
    solveEquationSet(&r0, &Margin_s[0], nDim, M);
    // r w u: reorder back
    r0.reorderElement(Index_invS);

    auto end = std::chrono::system_clock::now();  // Time end
    std::chrono::duration<double> elapsed_seconds = end - start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    Rcpp::Rcout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;


    /* ****************************************************************************************
    // Send output & free memory
    ***************************************************************************************** */
    // Output exp(r), exp(w), exp(u)
    for (k = 0, iDim = 0; iDim < nDim; iDim++) {
        for (i = 0; i < M[iDim]; i++) {
            R_out[k] = r0.r[iDim][i];
            k++;
        }
    }


    // Free memory
    for (iDim = 0; iDim < nDim; iDim++) {
        FREE(Margin[iDim]);
        FREE(Margin_s[iDim]);
        FREE(Index_invS[iDim]);
        FREE(Index_s[iDim]);
    }
    FREE(Y);
    FREE(Miss_);
    R_out.~NumericVector();

    // Return
    Result["ErrorFlag"] = 0;
    return Result;
}




void solveEquationSet(R_DIM* p_r, double** Margin, const uint nDim, const uint* M)
{
    // Parameters
    double errT_R = 1e-7;
    uint nTry_Max = 30;

    //
    uint iDim = 0, iTry = 0, curDim = 0;
    double diff = 100.0;

    R_DIM r_pre;
    r_pre = *p_r;
    r_pre.r[0][0] = p_r->r[0][0] + 100.0; // to make sure the initial diff is not smaller than the threshold




    //solveEquationSet_1Dim(p_r, 0, Margin, nDim, M);
    for (iTry = 0; iTry < nTry_Max; iTry++) {
        //iIter_G = iTry;  // DEBUG
        for (curDim = 0; curDim < nDim; curDim++) {
            //iDim_G = curDim; // DEBUG
            solveEquationSet_1Dim(iTry, p_r, curDim, Margin, nDim, M);
            //solveEquationSet_1Dim(iTry, p_r, 2, Margin, nDim, M);
        }
        diff = r_pre - *p_r;
        Rcpp::Rcout << "Iteration=" << iTry << ";  diff=" << diff << std::endl;
        if (diff < errT_R)
            break;
        r_pre = *p_r;
    }

}


//void solveEquationSet_1Dim(double* r, double** r_Part, const double* Margin_d, const uint nDim_Part, const uint* M_Part, const uint M_d)
void solveEquationSet_1Dim(int iIter_outer, R_DIM* p_r, const uint curDim, double** Margin, const uint nDim, const uint* M)
{
    // thld
    //double errT_r_base = 3.75e-10;
    double errT_r_base = 3.75e-9;

    //
    const double* Margin_d = Margin[curDim];
    const uint M_d = M[curDim];
    double* r_CurD = p_r->r[curDim];
    uint nDim_Part = nDim - 1;

    uint M_Part[MAX_DIM - 1] = { 0 };
    double* r_Part[MAX_DIM - 1];
    for (uint iDim = 0, kDim = 0; iDim < nDim; iDim++) {
        if (iDim == curDim) continue;
        r_Part[kDim] = p_r->r[iDim];
        M_Part[kDim] = M[iDim];
        kDim++;
    }

    // solve r
 //   for (uint i = 0; i < M_d; i++)
 //       r_CurD[i] = solve_r(i != 0, r_CurD[i], r_Part, nDim_Part, M_Part, Margin_d[i]);


    int i = 0, j = 0;
    double r_tmp = 0.0;
    double errT_r = 0.0;

    switch (iIter_outer) {
    case 0:
        errT_r = errT_r_base * 10000.0;
        break;
    case 1:
        errT_r = errT_r_base * 100.0;
        break;
    default:
        errT_r = errT_r_base;
    }


    int k = 0;
    int baseTmp = 1;
    for (i = 0; i < curDim; i++) {
        idx_Y_based[i] = baseTmp;
        baseTmp = baseTmp * M[i];
    }
    idx_Y_based[nDim] = baseTmp;
    baseTmp = baseTmp * M[curDim];
    for (i = curDim + 1; i < nDim; i++) {
        idx_Y_based[i - 1] = baseTmp;
        baseTmp = baseTmp * M[i];
    }



    if (iIter_outer < 2) { // First two iterations
        for (i = 0; i < M_d; i++) {
            idx_Y_based[nDim - 1] = idx_Y_based[nDim] * i;
            //if (iDebug == 1) printf("solveEquationSet_1Dim: direction = STARTING, curDim = %d, curElement = %d", curDim, i);
            if (iDebug == 1) Rcpp::Rcout << "solveEquationSet_1Dim: direction = STARTING, curDim = " << curDim << ", curElement = " << i << std::endl;
            r_tmp = solve_r(STARTING, i != 0, r_CurD[i], r_Part, nDim_Part, M_Part, Margin_d[i], M_d, errT_r);
            p_r->update(curDim, i, r_tmp);
        }
    }
    else { // Other iterations
        for (j = M_d - 1; j >= 0 && p_r->delta_r[curDim][j] >= 0; j--) {
            idx_Y_based[nDim - 1] = idx_Y_based[nDim] * j;
            //if (iDebug == 1) printf("solveEquationSet_1Dim: direction = ASCEND, curDim = %d, curElement = %d", curDim, j);
            if (iDebug == 1) Rcpp::Rcout << "solveEquationSet_1Dim: direction = ASCEND, curDim = " << curDim << ", curElement = " << j << std::endl;
            r_tmp = solve_r(ASCEND, j != M_d - 1, r_CurD[j], r_Part, nDim_Part, M_Part, Margin_d[j], M_d, errT_r);
            p_r->update(curDim, j, r_tmp);
        }
        for (i = j; i >= 0; i--) {
            //if (p_r->delta_r[curDim][i] > 0) 
            //    printf("p_r->delta_r[curDim][i] > 0!!! \n");
            idx_Y_based[nDim - 1] = idx_Y_based[nDim] * i;
            //if (iDebug == 1) printf("solveEquationSet_1Dim: direction = DESCEND, curDim = %d, curElement = %d", curDim, i);
            if (iDebug == 1) Rcpp::Rcout << "solveEquationSet_1Dim: direction = DESCEND, curDim = " << curDim << ", curElement = " << i << std::endl;
            r_tmp = solve_r(DESCEND, i != 0, r_CurD[i], r_Part, nDim_Part, M_Part, Margin_d[i], M_d, errT_r);
            p_r->update(curDim, i, r_tmp);
        }
    }

}






double solve_r(int Direction, bool iSameDim, double r_pre, double** r_Part, const uint nDim_Part, const uint* M_Part, const double m, const uint M_cur, double errT_r)
{
    // Parameters
    uint nTry_Max = 100;
    //double errT_r = 3.75e-10;
    uint thld_M_cur = 100; // If not sure which information to be borrowed, 
                            //then borrow from neighbor r_i if the length of current dimension > thld_M_cur
                            // otherwise, borrow from previous r_i(t-1)

    // Static variables (to stord the final results of last call)
    static double r_last = 0.0;
    static double f_last = 0.0;
    static double df_last = 0.0;
    static double m_last = 0.0;

    // Variables
    double r_try = 0.0;
    double f = 0.0, df = 0.0;
    double r_new = 0.0;
    uint iTry = 0;
    char iDir_inf = 0; // Indicate which information to be borrowed (0:neighbor r_i'; 1: previous r_i(t-1))



    // If the row or col are all 0 or 1, then directly assign r = 0 or r = Inf 
    if (m < 1e-8)
        return 0;
    double nEntries = 1;
    for (int i = 0; i < nDim_Part; i++)
        nEntries *= M_Part[i];
    if (abs(m - nEntries) < 1e-8 && (flag_Y_unique <= 2))  //!!! haven't consider the case where flag_Y_unique == 3; To do !!!
        return INF;



    // Decide which information to be borrowed
    r_try = r_last;
    f = f_last + m_last - m;
    df = df_last;
    r_new = r_try - f / df;
    r_new = r_new > 0 ? r_new : 0;

    if (iSameDim == 0)
        iDir_inf = 1;
    else {
        switch (Direction) {
        case ASCEND:
            if ((abs(r_new - r_try) > errT_r) && r_pre > r_new)
                iDir_inf = 1;
            break;
        case DESCEND:
            if (abs(r_new - r_try) > errT_r) {
                if (r_pre >= r_last || M_cur > thld_M_cur)
                    iDir_inf = 0;
                else if (abs(r_pre - r_new) < errT_r)
                    iDir_inf = 2; //Although abs(r_new - r_try) > errT_r, r_new is close enough
                else
                    iDir_inf = 1;
            }
            break;
        default:
            iDir_inf = 0;
        }
    }

    // Initialize the r_try & do first try
    if (flag_Miss == 1 || flag_Y_unique == 0) {
        if (iDir_inf == 1)
            r_try = r_pre;
        else
            r_try = r_new;
        f_df(&f, &df, r_try, r_Part, nDim_Part, M_Part, m);
        r_new = r_try - f / df;
        r_new = r_new > 0 ? r_new : 0;
    }
    else {
        if (iDir_inf == 1) {
            r_try = r_pre;
            f_df(&f, &df, r_try, r_Part, nDim_Part, M_Part, m);
            r_new = r_try - f / df;
            r_new = r_new > 0 ? r_new : 0;
        }
    }


    // Newton method to search the root:r_i
    if (abs(r_new - r_try) > errT_r && iDir_inf != 2) {
        r_try = r_new;
        for (iTry = 1; iTry < nTry_Max; iTry++) {
            f_df(&f, &df, r_try, r_Part, nDim_Part, M_Part, m);
            r_new = r_try - f / df;
            r_new = r_new > 0 ? r_new : 0;
            if (abs(r_new - r_try) < errT_r)
                break;
            r_try = r_new;

        }
    }



    if (iDebug == 1)
        printf("  nTry =%d,  r=%.8f\n", iTry, r_new);

    if (iTry >= nTry_Max - 1) {
        printf("Warning: when solving r, the number of iterations exceeds the maximun allowed iteration counts.\n");
    }

    // Stord the final results
    r_last = r_try;
    f_last = f;
    df_last = df;
    m_last = m;

    return r_new;
}






void f_df(double* f, double* df, double r_try, double** r_Part, const uint nDim_Part, const uint* M_Part, const double m)
{

    int iDim = 0;
    volatile char iExcess = 0;
    double A = 0, dA = 0;
    double prod_p = 1.0, prod = 1.0, prod_add1 = 1.0, inv_Prod_add1 = 1.0;
    uint y = 1;
    unsigned long long int idxTmp = 1;
    unsigned long long int nNonMiss = 0;
    int flag_INF = 0;

    //Index_HD idx_Part(nDim_Part, M_Part);
    //idx_Part.reset();

    int idx_Part[MAX_DIM] = { 0 };

    iExcess = 0;
    while (iExcess == 0) {
        flag_INF = 0;

        if (flag_Miss == 1 || flag_Y_unique == 0)
            for (idxTmp = idx_Y_based[nDim_Part], iDim = 0; iDim < nDim_Part; iDim++)
                idxTmp += idx_Y_based[iDim] * idx_Part[iDim];

        if (flag_Miss == 0 || Miss_[idxTmp]) {

            prod_p = 1.0;
            for (iDim = 0; iDim < nDim_Part; iDim++) {
                if (r_Part[iDim][idx_Part[iDim]] > INF_LB) {
                    flag_INF = 1;
                    break;
                }
                prod_p = prod_p * r_Part[iDim][idx_Part[iDim]];
            }


            prod = prod_p * r_try;
            prod_add1 = prod + 1;
            inv_Prod_add1 = 1.0 / prod_add1;


            if (flag_INF == 1) { // Cases where some of the w or u are Inf (due to the all 1 columns/pages)
                if (flag_Y_unique == 0) {
                    y = Y[idxTmp];
                    A += y;
                    dA += 0;
                }
                else {
                    A += 1;
                    dA += 0;
                }
            }
            else { // Normal cases
                if (flag_Y_unique == 0) {
                    y = Y[idxTmp];
                    A += y * prod * inv_Prod_add1;
                    dA += y * prod_p * inv_Prod_add1 * inv_Prod_add1;
                }
                else {
                    A += prod * inv_Prod_add1;
                    dA += prod_p * inv_Prod_add1 * inv_Prod_add1;
                }
            }
        }


        // idx_Part++;
        iExcess = 1;
        for (iDim = nDim_Part - 1; iDim >= 0; iDim--) {
            if (++idx_Part[iDim] < M_Part[iDim]) {
                iExcess = 0;
                break;
            }
            else {
                idx_Part[iDim] = 0;
            }
        }

    }

    *f = A - m;
    *df = dA;
}
