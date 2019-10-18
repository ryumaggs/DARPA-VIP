//mex -g -DDOUBLE covertree.C ThreadsWithCounter.C IDLList.C IDLListNode.C Vector.C Cover.C Point.C CoverNode.C EnlargeData.C Timer.C -lpthread

/*
 On OS/X, I had to run
 mex -setup
 from the MATLAB command line first. This creates the file
 ~/.matlab/mexopts.sh
 I then modified this file, in particular I set
 CC='gcc'
 in order to use gcc and not a specified version of gcc (as in the mexopts.sh file created by MATLAB) and
 SDKROOT='/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/'
 instead of the one created by MATLAB as things have moved with the latest version of Xcode (and Matlab also does not seem to properly recognized the latest version of OS X).
 */

/*
 outparams[0]=vectors.getIndex(cover.getRoot()->getPoint());
 outparams[1]=cover.getMinLevel();
 outparams[2]=cover.getNumLevels();
 outparams[3]=cover.getCount();
 outparams[4]=cover.getNumberInserted();
 outparams[5]=cover.getNumberDeep();
 outparams[6]=cover.getNumberDuplicates();
 outparams[7]=DISTANCE_FCN;
 */

#include "mex.h"
#include <math.h>
#include "Cover.H"
#include "CoverForest.H"
#include "string.h"
#include "Vector.H"
#include "ThreadsWithCounter.H"
#include "EnlargeData.H"
#include "Timer.H"
#include "VectorClassNames.H"
#include "covertree_MEXint.H"


//#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
//#define MAX(X, Y) (((X) < (Y)) ? (X) : (Y))

#define VERBOSE 10

int  dim=0;
REAL theta=0.0;
REAL mu=0.0;

mwSize dims[2];
int ndims=2;

const char* fnames_in[]={
    "theta",
    "numlevels",
    "minlevel",
    "NTHREADS",
    "BLOCKSIZE",
    "distancefcn",
    "classname",
    "NTREES"
};

const char* fnames_out_forest[]={
    "coverTrees",
    "ptIdxs",
    "nTrees"
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Distance_Mode DISTANCE_FCN = EUCLIDEAN;
    VectorClassNames VECTOR_CLASS;
    
#ifdef DOUBLE
    char mxIsClassName[] = "double";
#else
    char mxIsClassName [] = "single";
#endif
    
    /* Check for proper number of arguments. */
    if (nrhs != 2)                                                          { mexErrMsgTxt("Two inputs required."); }
    
    int nfields = mxGetNumberOfFields(prhs[0]);
    if(nfields<7)                                                           { mexErrMsgTxt("First input must have at least 7 fields."); }
    
    
    // First input: structure of options for the cover tree construction
    // Get theta
    const mxArray* tmp=mxGetField(prhs[0],0,fnames_in[0]);
    if(mxGetClassID(tmp)!=mxREAL_CLASS)                                     { mexErrMsgTxt("input.theta must be double\n"); }
    theta=*(REAL*)mxGetData(tmp);
    bool val=(theta>0.0)&&(theta<1.0);
    if(!val)                                                                { mexErrMsgTxt("Bad theta\n"); }
    mu=1.0/(1.0-theta);
    
    // Get outparams
    tmp=mxGetField(prhs[0],0,fnames_in[1]);
    if(mxGetClassID(tmp)!=mxINT32_CLASS)                                    { mexErrMsgTxt("input.numlevels must be int32\n"); }
    int numlevels=*(int*)mxGetData(tmp);
    
    // Get minlevel
    tmp=mxGetField(prhs[0],0,fnames_in[2]);
    if(mxGetClassID(tmp)!=mxINT32_CLASS)                                    { mexErrMsgTxt("input.minlevel must be int32\n"); }
    int minlevel=*(int*)mxGetData(tmp);
    
    // Get NTHREADS
    tmp=mxGetField(prhs[0],0,fnames_in[3]);
    if(mxGetClassID(tmp)!=mxINT32_CLASS)                                    { mexErrMsgTxt("input.NTHREADS must be int32\n"); }
    int NTHREADS=*(int*)mxGetData(tmp);
    
    // Get BLOCKSIZE
    tmp=mxGetField(prhs[0],0,fnames_in[4]);
    if(mxGetClassID(tmp)!=mxINT32_CLASS)                                    { mexErrMsgTxt("input.BLOCKSIZE must be int32\n"); }
    int BLOCKSIZE=*(int*)mxGetData(tmp);
    
    // Get distancefcn
    tmp=mxGetField(prhs[0],0,fnames_in[5]);
    if(tmp==NULL)   {
        DISTANCE_FCN = EUCLIDEAN;
    }
    else {
        if(mxGetClassID(tmp)!=mxINT32_CLASS)                                { mexErrMsgTxt("input.distancefcn must be int32\n"); }
        else    {
            int *pDISTANCE_FCN = (int*)mxGetData(tmp);
            DISTANCE_FCN = (Distance_Mode)(*pDISTANCE_FCN);
        }
    }
    // Get classname
    tmp=mxGetField(prhs[0],0,fnames_in[6]);
    if(tmp==NULL)   {
        VECTOR_CLASS = VECTORS;
    }
    else {
        if(mxGetClassID(tmp)!=mxINT32_CLASS)                                { mexErrMsgTxt("input.classname must be int32\n"); }
        VECTOR_CLASS = *(VectorClassNames*)mxGetData(tmp);
    }
    
    // Get NTHREADS
    tmp=mxGetField(prhs[0],0,fnames_in[7]);
    if(mxGetClassID(tmp)!=mxINT32_CLASS)                                     { mexErrMsgTxt("input.NTREES must be int32\n"); }
    int* pNTREES=(int*)mxGetData(tmp);
    int NTREES=(int)*pNTREES;
    
    
    // Second Input: X vectors
    Vectors *vectors;
    LoadXToVectors( prhs[1], &vectors, VECTOR_CLASS, DISTANCE_FCN );
    int N   = vectors->getCount();
    int dim = vectors->getDim();
    
    // Create covertree
    if( N==1 ) {    NTHREADS = 0;   }
    
    ThreadsWithCounter threads(NTHREADS);
    
    if( NTREES<=1 ) {
        if( VERBOSE>=1 )                                                    { mexPrintf("\nConstructing cover tree..."); }
        // Construct covertree
        SegList< DLPtrListNode<CoverNode> > seglist(N);
        const       Vector* vector=(vectors->next());
        EnlargeData enlargedata(&threads,BLOCKSIZE,vectors->getRemaining());
        
        Cover       cover(vector,seglist,numlevels,minlevel);
        if( VERBOSE>=2 )                                                    { mexPrintf("\nEnlarging cover tree..."); }
        cover.enlargeBy(enlargedata,*vectors);
        
        // Create the return argument
        if( VERBOSE>=2 )                                                    { mexPrintf("\nConstructing return argument..."); }
        SaveCoverTreeToMatStruct( plhs, *vectors, cover, &enlargedata, VECTOR_CLASS );
        if( VERBOSE>=1 )                                                    { mexPrintf("done"); }
    } else {
        if( VERBOSE>=1 )                                                    { mexPrintf("\nConstructing cover forest..."); }
        // Construct cover forest
        EnlargeData **enlargedata_forest;
        CoverForest coverForest( NTREES );
        if( VERBOSE>=2 )                                                    { mexPrintf("\nEnlarging cover forest..."); }
        coverForest.enlargeBy(*vectors, N, dim, numlevels, BLOCKSIZE, &enlargedata_forest );
        
        // Create the return argument
        if( VERBOSE>=2 )                                                    { mexPrintf("\nConstructing return argument..."); }
        mxArray* plhs_covertrees=0;
        for( INDEX i=0; i<NTREES; i++ ) {
            SaveCoverTreeToMatStruct( &plhs_covertrees, *vectors, *(coverForest.getCoverTree(i)), enlargedata_forest[i], VECTOR_CLASS, NTREES, i, coverForest.getptIdxs()[2*i]);
        }
        *plhs = mxCreateStructMatrix(1, 1, 3, fnames_out_forest);
        mxSetField(plhs[0],0,fnames_out_forest[0], plhs_covertrees );
        
        if( VERBOSE>=2 )                                                    { mexPrintf("\tplhs.%s...",fnames_out_forest[1]); }
        ndims           = 2;                                                                                                                                                                    // ptIdxs
        dims[0]         = 2;    dims[1]         = NTREES;
        mxArray *fout   = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
        mxSetField(plhs[0],0,fnames_out_forest[1], fout);
        memcpy(mxGetData(fout),coverForest.getptIdxs(),2*NTREES*sizeof(unsigned int));
        
        if( VERBOSE>=2 )                                                    { mexPrintf("\tplhs.%s...",fnames_out_forest[2]); }
        ndims           = 2;                                                                                                                                                                    // NTREES
        dims[0]         = 1;    dims[1]         = 1;
        fout            = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
        *(int*)mxGetData(fout)   = NTREES;
        mxSetField(plhs[0],0,fnames_out_forest[2], fout);
        if( VERBOSE>=1 )                                                    { mexPrintf("done"); }
    }
    
    // Clean up
    delete vectors;
}


#include "FastSeg.C"
template class SegList<DLPtrListNode<CoverNode> >;
