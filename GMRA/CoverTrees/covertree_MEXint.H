//
//  covertree_MEXint.hpp
//  
//
//  Created by Mauro Maggioni on 1/26/17.
//
//

#ifndef covertree_MEXint_hpp
#define covertree_MEXint_hpp

#include <stdio.h>
#include "mex.h"
#include <math.h>
#include "Cover.H"
#include "string.h"
#include "VectorClassNames.H"
#include "Vector.H"


#ifdef DOUBLE
#define mxREAL_CLASS mxDOUBLE_CLASS
#else
#define mxREAL_CLASS mxSINGLE_CLASS
#endif


extern const char* covertreestruct_names[];
extern REAL theta;
extern REAL mu;

mxArray* SaveCoverTreeToMatStruct( mxArray **plhs, const Vectors& vectors, const Cover& cover, const EnlargeData *enlargedata, VectorClassNames VECTOR_CLASS, int structArrayLen=1, mwIndex idx_plhs=0, unsigned int indexoffset = 0);

int LoadXToVectors( const mxArray *mX, Vectors** vectorsX, VectorClassNames VECTOR_CLASS, Distance_Mode DISTANCE_FCN = EUCLIDEAN );
int LoadMatStructToCoverTree( const mxArray *prhs[], Distance_Mode* DISTANCE_FCN, VectorClassNames* VECTOR_CLASS, int* minlevel, int* numlevels, int **levels );

#endif /* covertree_MEXint_hpp */
