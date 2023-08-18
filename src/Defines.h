#pragma once


#ifndef Defines_H
#define Defines_H



/************* HyperParameter define ******************/
#define INF 1e+100
#define INF_LB 1e+99 // any value larger than this value will be regarded as INF
#define MAX_COL_CSV 10
#define MAX_ROW_CSV 30000

#define MAX_DIM 6 // max dimensions



#define ASCEND 1
#define DESCEND 0
#define NODIRECTION 2
#define STARTING 3

/************* Type define ******************/
typedef unsigned int uint;







/********** Macros for memory menagement **********/

#define MALLOC(type) (type *)malloc((sizeof(type)))
#define CALLOC(type,number) (type *)calloc((number),(sizeof(type)))
#define FREE(pointer) { if((pointer)!=NULL) free((void *)pointer); }


#endif