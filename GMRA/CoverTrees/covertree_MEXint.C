//
//  covertree_MEXint.cpp
//
//
//  Created by Mauro Maggioni on 1/26/17.
//
//

#include "covertree_MEXint.H"
#include "mex.h"
#include <math.h>
#include "CoverForest.H"
#include "string.h"
#include "ThreadsWithCounter.H"
#include "EnlargeData.H"
#include "Timer.H"

#define DEBUG

const char* covertreestruct_names[]={
    "theta",
    "outparams",
    "radii",
    "levels",
    "ncallstogetdist"
};


// Converts covertree to MATLAB structure
// This is designed for saving a single covertree, but the parameters idx_plhs and indexoffset are designed to help saving a coverforest, by saving multiple covertrees in a struct array. idx_plhs
// is the index into the struct array (typically equal to the index of the covertree in the coverforest) and indexoffset is used when creating CoverIndices for the covertrees in a coverforest, since
// the i-th covertree in a coverforest is about points with indices given by coverForest.ptIdxs(:,i). indexoffset adjusts for this.
mxArray* SaveCoverTreeToMatStruct( mxArray **plhs, const Vectors& vectors, const Cover& cover, const EnlargeData *enlargedata, VectorClassNames VECTOR_CLASS, int structArrayLen, mwIndex idx_plhs, int indexoffset)   {
    
    if( *plhs==0 )  { *plhs = mxCreateStructMatrix(1, structArrayLen, 5, covertreestruct_names); }
    
    mxArray* fout;
    mwSize ndims = 2;
    int dims[ndims];
    
    dims[0]=1; dims[1]=1;                                                                                                                                                                       // theta
    fout    = mxCreateNumericArray(ndims,dims,mxREAL_CLASS,mxREAL);
    REAL* p = (REAL*)mxGetData(fout);
    mxSetField(*plhs,idx_plhs,covertreestruct_names[0],fout);
    p[0]    = theta;
    
    dims[0]=1; dims[1]=9;                                                                                                                                                                       // outparams
    fout    = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int* outparams=(int*)mxGetData(fout);
    mxSetField(*plhs,idx_plhs,covertreestruct_names[1], fout);
    outparams[0] = vectors.getIndex(cover.getRoot()->getPoint());
    outparams[1] = cover.getMinLevel();
    outparams[2] = cover.getNumLevels();
    outparams[3] = cover.getCount();
    outparams[4] = cover.getNumberInserted();
    outparams[5] = cover.getNumberDeep();
    outparams[6] = cover.getNumberDuplicates();
    outparams[7] = cover.getDistanceMode();
    outparams[8] = VECTOR_CLASS;
    
    dims[0]=1; dims[1]=cover.getNumLevels();                                                                                                                                                    // radii
    fout         = mxCreateNumericArray(ndims,dims,mxREAL_CLASS,mxREAL);
    REAL* pradii =(REAL*)mxGetData(fout);
    pradii[0]    = cover.getMaxRadius();
    for(int i=1;i<cover.getNumLevels();i++) {   pradii[i]=theta*pradii[i-1];    }
    mxSetField(*plhs,idx_plhs,covertreestruct_names[2], fout);
    
    dims[0]     = cover.getCount(); dims[1]=5;                                                                                                                                                      // levels
    mexPrintf("\n cover(%d).getCount()=%d",idx_plhs,dims[0]);
    mexPrintf("\n indexoffset=%d",indexoffset);
    fout        = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    int* base   = (int*)mxGetData(fout);
    CoverIndices coverindices(&cover,&vectors,base,indexoffset);
#ifdef DEBUG
    mexPrintf("\n cover(%d)=%p",idx_plhs,&cover);
    mexPrintf("\n cover(%d)->first()=%p",idx_plhs,cover.first());
    mexPrintf("\n cover(%d)->first()->getPoint()->getIndex=%d",idx_plhs,vectors.getIndex(cover.first()->getPtr()->getPoint()));
    mexPrintf("\n cover(%d)->first()->getLevel()=%d",idx_plhs,cover.first()->getPtr()->getLevel());
    mexPrintf("\n cover(%d)->first()->next()->getPoint()->getIndex=%d",idx_plhs,vectors.getIndex(cover.first()->next()->getPtr()->getPoint()));
    mexPrintf("\n cover(%d)->first()->next()->getLevel()=%d",idx_plhs,cover.first()->next()->getPtr()->getLevel());
    mexPrintf("\n base:\n");
    for( unsigned int i = 0; i<4; i++) {
        for( unsigned int j = 0; j<5; j++ ) {
            mexPrintf("%d ",base[i*5+j]);
        }
    }
#endif
    mxSetField(*plhs,idx_plhs,covertreestruct_names[3], fout);
    
    dims[0]=1; dims[1]=2;                                                                                                                                                                       // ncallstogetdist
    fout = mxCreateNumericArray(ndims,dims,mxINT64_CLASS,mxREAL);
    long int* pncalls=(long int*)mxGetData(fout);
    mxSetField(*plhs,idx_plhs,covertreestruct_names[4], fout);
    pncalls[0]  = enlargedata[0].getMergeNCallsToGetDist();
    pncalls[1]  = enlargedata[0].getThreadNCallsToGetDist();
    
    return 0;
}


int LoadXToVectors( const mxArray *mX, Vectors** vectorsX, VectorClassNames VECTOR_CLASS, Distance_Mode DISTANCE_FCN )   {
    mwSize ndims_in=mxGetNumberOfDimensions(mX);
    const mwSize* dims_in=mxGetDimensions(mX);
    if(mxGetClassID(mX)!=mxREAL_CLASS)                                                          { mexErrMsgTxt("Second input is the wrong type float/double\n"); return -1; }
    
    int N = dims_in[ndims_in-1];                                                                                                                                                                // Number of points
    int dim   = 1;                                                                                                                                                                              // Ambient dimension
    for(int i=0;i<ndims_in-1;i++) { dim*=dims_in[i]; }
    
    // Create Vectors
    REAL* X = (REAL*)mxGetData(mX);
    switch (VECTOR_CLASS)   {
        case VECTORS:
            *vectorsX = new Vectors(X,N,dim,DISTANCE_FCN);
            break;
        case IMAGES:
            *vectorsX = new Images(X,N,sqrt(dim),sqrt(dim),DISTANCE_FCN);
            break;
        case MOLECULARSTATES:
            if( dim % 3==0 )                                                                    {   *vectorsX = new MolecularStates(X,N,dim/3);     }
            else                                                                                {   mexErrMsgTxt("Number of dimensions should be a multiple of 3\n");  return -1; }
            break;
        default:
            mexErrMsgTxt("\n Invalid classname for Vectors.");
            break;
    }
    
    return 0;
}




// Converts Matlab structure to covertree
int LoadMatStructToCoverTree( const mxArray *prhs[], Distance_Mode* DISTANCE_FCN, VectorClassNames* VECTOR_CLASS, Vectors **vectorsX, Cover **cover, SegList< DLPtrListNode<CoverNode> > **seglist, int indexoffset ) {
    int minlevel, numlevels, *levels;

#ifdef DOUBLE
    char mxIsClassName[] = "double";
#else
    char mxIsClassName [] = "single";
#endif

    // Check the first argument, which is what was returned by covertree
    int nelements=mxGetNumberOfFields(prhs[0]);
    if(nelements!=5)                                                                        { mexErrMsgTxt("Covertree structure should have five fields.");     return -1; }
    
    const mxArray* tmp=mxGetField(prhs[0],0,covertreestruct_names[0]);                                                                                                                          // theta
    if(!mxIsClass(tmp,mxIsClassName))                                                       { mexErrMsgTxt("First field of Covertree structure must be of correct real type\n"); return -1; }
    REAL* ptheta=(REAL*)mxGetData(tmp);
    theta=*ptheta;
    mu=1.0/(1.0-theta);
    
    // Get second field of first input
    tmp=mxGetField(prhs[0],0,covertreestruct_names[1]);                                                                                                                                         // outparams
    if(!mxIsClass(tmp,"int32"))                                                             { mexErrMsgTxt("Second field of Covertree structure must be int32\n"); return -1; }
    // Get dimensions of second field of first input
    mwSize ndims_in=mxGetNumberOfDimensions(tmp);
    const mwSize* dims_in=mxGetDimensions(tmp);
    bool val=(ndims_in==2)&&(dims_in[0]==1)&&(dims_in[1]==9);
    if(!val)                                                                                { mexErrMsgTxt("Second field of Covertree structure has bad size\n"); return -1; }
    int* outparams  = (int*)mxGetData(tmp);
    minlevel    = outparams[1];
    numlevels   = outparams[2];
    *DISTANCE_FCN   = (Distance_Mode)(outparams[7]);
    *VECTOR_CLASS   = (VectorClassNames)(outparams[8]);
    
    tmp=mxGetField(prhs[0],0,covertreestruct_names[3]);                                                                                                                                         // levels
    if(!mxIsClass(tmp,"int32"))                                                             { mexErrMsgTxt("Fourth field of Covertree structure should be int32\n"); return -1; }
    mwSize ndims_indices = mxGetNumberOfDimensions(tmp);
    const mwSize* dims_indices = mxGetDimensions(tmp);
    if(ndims_indices!=2)                                                                    { mexErrMsgTxt("Fourth field of Covertree structure has bad size\n");   return -1; }
    levels=(int*)mxGetData(tmp);
    
    // Process the second field of input, which are the vectors X
    tmp=prhs[1];
    if(!mxIsClass(tmp,mxIsClassName))                                                       { mexErrMsgTxt("Second input must be of correct real type\n"); }
    mwSize ndims_X=mxGetNumberOfDimensions(tmp);
    const mwSize* dims_X=mxGetDimensions(tmp);
    

//Replace with: int LoadXToVectors( const mxArray *mX, Vectors** vectorsX, VectorClassNames VECTOR_CLASS, Distance_Mode DISTANCE_FCN )
    
    
    int NX=dims_X[ndims_in-1];
    int d=dims_indices[0];
    if ( d!=NX )                                                                            { mexErrMsgTxt("Mismatch between first and second inputs\n"); }
    int dim=1;
    for(int i=0;i<ndims_X-1;i++) {
        dim*=dims_X[i];
    }
    
    // Construct vectors X
    REAL* X=(REAL*)mxGetData(tmp);
    switch (*VECTOR_CLASS)   {
        case VECTORS:
            *vectorsX = new Vectors(X,NX,dim,DISTANCE_FCN);
            break;
        case IMAGES:
            *vectorsX = new Images(X,NX,sqrt(dim),sqrt(dim),DISTANCE_FCN);
            break;
        case MOLECULARSTATES:
            if( dim % 3==0 )    {
                *vectorsX = new MolecularStates(X,NX,dim/3);
            }   else    {
                mexErrMsgTxt("Number of dimensions should be a multiple of 3\n");
            }
            break;
        default:
            mexErrMsgTxt("\n Invalid classname.");
            break;
    }
    
    // Re-format covertree
    CoverIndices coverindices(theta,numlevels,minlevel,NX,indices,-indexoffset);
    *seglist = new SegList< DLPtrListNode<CoverNode> >(16384);
    *cover = new Cover(**vectorsX,*seglist,coverindices);
    
    return 0;
}
