// CMC_V0.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <math.h>       /* exp */
#include <string.h>

#include "CMC_V0.h"

#include "r_solver.cpp"
// #include "ReadCSV.h"

// running time
#include <chrono>
#include <ctime>  

#include <algorithm> // for sort
//#include <vector>  // for vector


#include "mex.h"
#include <matrix.h>



void solveEquationSet(R_DIM* p_r, double** Margin, const uint nDim, const uint* M);
void solveEquationSet_1Dim(int iIter_outer, R_DIM* p_r, const uint curDim, double** Margin, const uint nDim, const uint* M);
double solve_r(int Direction, bool iSameDim, double r_pre, double** r_Part, const uint nDim_Part, const uint* M_Part, const double m, const uint M_cur, double errT_r = 3.75e-10);
void f_df(double* f, double* df, double r_try, double** r_Part, const uint nDim_Part, const uint* M_Part, const double m);



unsigned int nIterA0 = 0;
unsigned int nIterA = 0;
unsigned int nIterA_Dim0 = 0;
unsigned int nIterA_Dim1 = 0;
unsigned int nIterA_Dim2 = 0;
unsigned int iDim_G = 0;
unsigned int iIter_G = 0;


uint* Y = NULL;
uint flag_Y_unique = 0;
unsigned long long int idx_Y_based[MAX_DIM + 1] = { 0 };
bool* Miss_ = NULL;
uint flag_Miss = 0;


void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[]) {

    int nInput = nrhs;
    int nOutput = nlhs;


    uint nDim = 0;
    double* Margin[MAX_DIM];  // Marginal total counts
    uint M[MAX_DIM] = { 0 }; //Length of each dimension

    uint len_Y = 0;
    //uint* Y; // global variable
    unsigned long long int len_Miss_ = 0;
    //bool* Miss_; // global variable
    //uint flag_Miss = 0;  // global variable

    int j = 0;

    // Get margin inf
    const mxArray* Margin_Cell;
    Margin_Cell = prhs[0];
    if (!mxIsCell(Margin_Cell))
        mexErrMsgTxt("Margin must be a cell.");

    nDim = mxGetNumberOfElements(Margin_Cell);
    printf("nDim =%d\n", nDim);
    int iDim = 0;
    mxArray* cellElement;
    for (iDim = 0; iDim < nDim; iDim++) {
        cellElement = mxGetCell(Margin_Cell, iDim);
        M[iDim] = mxGetNumberOfElements(cellElement);
        Margin[iDim] = mxGetPr(cellElement);
        printf("Dim=%d has %d elements.\n", iDim, M[iDim]);
    }


    // Get Y
    // Y: 1 X 1 or  N X 1; where N is the total elements of the tensor
    if (nInput > 1) {
        const mxArray* Y_Vec;
        Y_Vec = prhs[1];
        len_Y = mxGetNumberOfElements(Y_Vec);
        if (len_Y == 0) {
            len_Y = 1;
            Y = MALLOC(uint);
            Y[0] = 1;
        }
        else
            Y = (uint*)mxGetPr(Y_Vec);
    }
    else {
        len_Y = 1;
        Y = MALLOC(uint);
        Y[0] = 1;
    }


    // Get Miss_ (i.e., indicate whether the element is NOT missing)
    // Miss_ : N X 1; where N is the total elements of the tensor
    
    if (nInput > 2) {
        const mxArray* Miss__Vec;
        Miss__Vec = prhs[2];
        len_Miss_ = mxGetNumberOfElements(Miss__Vec);
        if (len_Miss_ == 0) {
            flag_Miss = 0;
        }
        else {
            flag_Miss = 1;
            Miss_ = (bool*)mxGetPr(Miss__Vec);
        }
    }
    else {
        flag_Miss = 0;
    }


    
    /*
    //printf("len_Miss_: %.1f.\n", (double)len_Miss_);
    if (flag_Miss == 1) {
        printf("Miss_ is:\n");
        for (j = 0; j < (len_Miss_ < 10? len_Miss_ : 10); j++)
            printf("%d ", Miss_[j]);
        printf("\n");

        printf("Miss_ is (end):\n");
        for (long long int jj = len_Miss_-10; jj < len_Miss_; jj++)
            printf("%d ", Miss_[jj]);
        printf("\n");

    }
    */


    //mwSize numDataPoints;
    //numDataPoints = mxGetNumberOfElements(prhs[2]);
    //int64_t tmp1 = numDataPoints;


    // Call CMC
    R_DIM r_out(nDim, M);
    CMC(&r_out,nDim, Margin,M,len_Y,len_Miss_);




    // return r
    double** result = r_out.r;
    double* retPtr;
    plhs[0] = mxCreateCellMatrix(nDim, 1);
    for (iDim = 0; iDim < nDim; iDim++) {
        mxArray* ithCell = mxCreateDoubleMatrix(M[iDim], 1, mxREAL);
        retPtr = mxGetPr(ithCell);
        for (j = 0; j < M[iDim]; j++)
            retPtr[j] = result[iDim][j];  // return values
        mxSetCell(plhs[0], iDim, ithCell);
    }

}


//int main()
//{
//    std::cout << "Hello World!2\n";
//    main_tmp();
//}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file



void CMC(R_DIM* r_out, uint nDim, double** Margin, uint* M, uint len_Y, unsigned long long int len_Miss_)
{
    //uint nDim = 2;
    uint iDim = 0, i = 0, j = 0, k = 0;

    /* Read marginal information */
    //char fileName_m_base[] = "L:/CMC_C_R_code/data/m0.bin";

    //double* Margin[MAX_DIM];  // Marginal total counts
    //uint M[MAX_DIM]= {0}; //Length of each dimension
    //size_t len = strlen(fileName_m_base);
    //for (iDim = 0; iDim < nDim; iDim++) {        
    //    readBin(&Margin[iDim], fileName_m_base, &M[iDim]);
    //   fileName_m_base[len-5] ++;
    //}


    // Check. The sum of margins in each dimension should be equal
    double sumMargin[MAX_DIM] = {0};
    for (iDim = 0; iDim < nDim; iDim++) {
        for (j = 0; j < M[iDim]; j++) {
            sumMargin[iDim] += Margin[iDim][j];
        }
    }
    for (i = 0; i < nDim; i++) {
        if (abs((sumMargin[i] - sumMargin[0])/sumMargin[0] ) > 0.0001) {
            std::cout << "ERROR! The sum of margins in each dimension should be equal.";
            exit(0);
        }
    }
    
    // Output the dimension information
    std::cout << "Number of dimensions=" << nDim << "\n";
    for (iDim = 0; iDim < nDim; iDim++)
        std::cout << M[iDim] << " ";
    std::cout << "\n";


    // Sorting marginal totals (for speed up in latter step)
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
        sort(Index_s[iDim], Index_s[iDim]+M[iDim], \
            [&](int ii, int jj) {return Margin_s[iDim][ii] < Margin_s[iDim][jj]; });
        sort(Index_invS[iDim], Index_invS[iDim] + M[iDim], \
            [&](int ii, int jj) {return Index_s[iDim][ii] < Index_s[iDim][jj]; });
        sort(Margin_s[iDim], Margin_s[iDim] + M[iDim]);
    }




    // Read Y
    uint Y0 = 1;
    unsigned long long int tmp_len = 1;
    
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
        // Check the length of Y
        for (iDim = 0; iDim < nDim; iDim++)
            tmp_len = tmp_len * M[iDim];
        if (tmp_len != len_Y) {
            std::cout << "ERROR! The length of Y is not consistent with the inputted margin.";
            exit(0);
        }
    }


    // Divide the margins by Y_ij, if Y_ij are all the same but not equal to 1
    if (flag_Y_unique == 1) {
        for (uint iDim = 0; iDim < nDim; iDim++) {
            for (j = 0; j < M[iDim]; j++) {
                Margin_s[iDim][j] = Margin_s[iDim][j] / Y0;
            }
        }
    }



     //Read Missing inf
    if (flag_Miss) {
        printf("Reading index of missing entries.\n");
        for (tmp_len = 1, iDim = 0; iDim < nDim; iDim++)
            tmp_len = tmp_len * M[iDim];

        if (tmp_len != len_Miss_) {
            printf( "ERROR! The length of Miss_ is not consistent with the inputted margin.");
            exit(0);
        }
    }
    else {
        printf("No entries are marked as missing.\n");
    }

    


    /* ***************************
    // Find r w u
    *************************** */
    nIterA = 0;
    nIterA0 = 0;
    nIterA_Dim0 = 0;
    nIterA_Dim1 = 0;
    nIterA_Dim2 = 0;

    
    R_DIM r0(nDim, M);

    //solveEquationSet;
    auto start = std::chrono::system_clock::now();
    solveEquationSet( &r0, &Margin_s[0], nDim, M);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    //std::cout <<  "elapsed time: " << elapsed_seconds.count() << "s\n";
    printf ( "elapsed time: %.5f s\n", elapsed_seconds.count());


    // r w u: reorder back
    *r_out = r0;
    (*r_out).reorderElement(Index_invS);


 

    // Free memory
    for (iDim = 0; iDim < nDim; iDim++) {
        //FREE(Margin[iDim]);
        FREE(Margin_s[iDim]);
        FREE(Index_invS[iDim]);
        FREE(Index_s[iDim]);
    }
    
    return;

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
        iIter_G = iTry;  // DEBUG
        for (curDim = 0; curDim < nDim; curDim++) {
            iDim_G = curDim; // DEBUG
            solveEquationSet_1Dim(iTry, p_r, curDim, Margin, nDim, M);
            //solveEquationSet_1Dim(iTry, p_r, 2, Margin, nDim, M);
        }
        diff = r_pre - *p_r;
        std::cout << "iTry=" << iTry << ";  diff=" << diff << "\n";
        if (diff < errT_R)
            break;
        r_pre = *p_r;
    }

}


//void solveEquationSet_1Dim(double* r, double** r_Part, const double* Margin_d, const uint nDim_Part, const uint* M_Part, const uint M_d)
void solveEquationSet_1Dim(int iIter_outer, R_DIM* p_r, const uint curDim, double** Margin, const uint nDim, const uint* M)
{
    // thld
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
            r_tmp = solve_r(STARTING, i != 0, r_CurD[i], r_Part, nDim_Part, M_Part, Margin_d[i], M_d, errT_r);
            p_r->update(curDim, i, r_tmp);
        }
    }
    else { // Other iterations
        for (j = M_d - 1; j >= 0 && p_r->delta_r[curDim][j] >= 0; j--) {
            idx_Y_based[nDim - 1] = idx_Y_based[nDim] * j;
            r_tmp = solve_r(ASCEND, j != M_d - 1, r_CurD[j], r_Part, nDim_Part, M_Part, Margin_d[j], M_d, errT_r);
            p_r->update(curDim, j, r_tmp);
        }
        for (i = j; i >= 0; i--) {
            //if (p_r->delta_r[curDim][i] > 0) 
            //    printf("p_r->delta_r[curDim][i] > 0!!! \n");
            idx_Y_based[nDim - 1] = idx_Y_based[nDim] * i;
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
    bool iDir_inf = 0; // Indicate which information to be borrowed (0:neighbor r_i'; 1: previous r_i(t-1))


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
    if (flag_Y_unique != 0) {
        if (iDir_inf == 1) {
            r_try = r_pre;
            f_df(&f, &df, r_try, r_Part, nDim_Part, M_Part, m);
            r_new = r_try - f / df;
            r_new = r_new > 0 ? r_new : 0;
        }
    }
    else {
        if (iDir_inf == 1)
            r_try = r_pre;
        else
            r_try = r_new;
        f_df(&f, &df, r_try, r_Part, nDim_Part, M_Part, m);
        r_new = r_try - f / df;
        r_new = r_new > 0 ? r_new : 0;
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




    nIterA0++;

    if (iTry >= nTry_Max - 1) {
        std::cout << "ERROR when solving r. The number of iterations exceeds the maximun allowed iteration counts.";
        exit(0);
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
    register double A = 0, dA = 0;
    register double prod_p = 1.0, prod = 1.0, prod_add1 = 1.0, inv_Prod_add1 = 1.0;
    int y = 1;
    int idxTmp = 1;


    //Index_HD idx_Part(nDim_Part, M_Part);
    //idx_Part.reset();

    register int idx_Part[MAX_DIM] = { 0 };

    iExcess = 0;
    while (iExcess == 0) {
        prod_p = 1.0;
        for (iDim = 0; iDim < nDim_Part; iDim++) {
            prod_p = prod_p * r_Part[iDim][idx_Part[iDim]];
        }
        prod = prod_p * r_try;
        prod_add1 = prod + 1;
        inv_Prod_add1 = 1.0 / prod_add1;

        if (flag_Y_unique == 0) {
            for (idxTmp = idx_Y_based[nDim_Part], iDim = 0; iDim < nDim_Part; iDim++)
                idxTmp += idx_Y_based[iDim] * idx_Part[iDim];
            y = Y[idxTmp];
            A += y * prod * inv_Prod_add1;
            dA += y * prod_p * inv_Prod_add1 * inv_Prod_add1;
        }
        else {
            A += prod * inv_Prod_add1;
            dA += prod_p * inv_Prod_add1 * inv_Prod_add1;
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

    nIterA++;

    switch (iIter_G) {
    case 0:
        nIterA_Dim0++;
        break;
    case 1:
        nIterA_Dim1++;
        break;
    case 2:
        nIterA_Dim2++;
        break;
    default:;
    }



}
