//Compile command: mex -O -DDOUBLE findnearest.C ThreadsWithCounter.C IDLList.C IDLListNode.C Vector.C Cover.C Point.C CoverNode.C EnlargeData.C FindWithinData.C findNearest.C -lpthread

/* findnearest take five arguments:
 The first argument is what was returned by covertree
 The second argument is the array that was passed as the second argument to
 covertree
 The third argument is the array of query vectors
 The fourth argument the number of nearest neighbors to find
 The fifth argument, of class int32, is the number of threads (0 for serial).
 
 findnearest returns an k x M matrix of class int32 whose j-th column
 is the indices of k nearest nearest neighbors of the j-th queryvector
 */

#include <math.h>
#include "mex.h"
#include "Cover.H"
#include "Vector.H"
#include "ThreadsWithCounter.H"
#include "FindNearestData.H"
#include "Distances.H"
#include "VectorClassNames.H"

const char* fnames_in[]={
    "theta",
    "outparams",
    "radii",
    "levels",
    "ncallstogetdist",
};

const char* fnames_out[]={
    "indices",
    "distances"
};

int dim=0;
REAL theta, mu;
int ndims=2;
mwSize dims[2];









void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Distance_Mode DISTANCE_FCN;
    VectorClassNames VECTOR_CLASS;
#ifdef DOUBLE
    char mxIsClassName[] = "double";
#else
    char mxIsClassName [] = "single";
#endif
    SegList< DLPtrListNode<CoverNode> > *seglist;
    if ( LoadMatStructToCoverTree( prhs[], &DISTANCE_FCN, &VECTOR_CLASS, &vectorsX, &cover, &seglist )!= 0 ) { mexErrMsgTxt("Invalid covertree data"); }
    
    // Construct vectors Y
    tmp=prhs[2];
    mwSize ndims_Y=mxGetNumberOfDimensions(tmp);
    const mwSize* dims_Y=mxGetDimensions(tmp);
    if(!mxIsClass(tmp,mxIsClassName))                        { mexErrMsgTxt("Third argument should be double\n"); }
    int NY=dims_Y[ndims_Y-1];
    int dimY=1;
    for(int i=0;i<ndims_Y-1;i++) {
        dimY*=dims_Y[i];
    }
    if(dimY!=dim)                                       { mexErrMsgTxt("Dimension mismatch\n"); }
    REAL* Y=(REAL*)mxGetData(tmp);
    Vectors *vectorsY;
    switch (VECTOR_CLASS)   {
        case VECTORS:
            vectorsY = new Vectors(Y,NY,dim,DISTANCE_FCN);            
            break;
        case IMAGES:
            vectorsY = new Images(Y,NY,sqrt(dim),sqrt(dim),DISTANCE_FCN);
            break;
        case MOLECULARSTATES:
            if( dim % 3==0 )    {
                vectorsY = new MolecularStates(Y,NY,dim/3);
            }   else    {
                mexErrMsgTxt("Number of dimensions should be a multiple of 3\n");
            }
            break;
        default:
            mexErrMsgTxt("\n Invalid classname.");
            break;
    }
    
    // Get number of nearest neighbors to compute
    tmp=prhs[3];
    ndims_in    = mxGetNumberOfDimensions(tmp);
    dims_in     = mxGetDimensions(tmp);
    if(!mxIsClass(tmp,"int32"))                         { mexErrMsgTxt("Fifth argument should be int32\n"); }
    int* pk=(int*)mxGetData(tmp);
    int k=*pk;
    
    // Get number of threads and create them 
    tmp=prhs[4];
    ndims_in    = mxGetNumberOfDimensions(tmp);
    dims_in     = mxGetDimensions(tmp);
    if(!mxIsClass(tmp,"int32"))                         { mexErrMsgTxt("Fifth argument should be int32\n"); }
    int NTHREADS= *(int*)mxGetData(tmp);
    
    //const Point** ptarr=(const Point**)mxMalloc((k+1)*NY*sizeof(const Point*));   double* distances=(double*)mxMalloc(k*NY*sizeof(double));
    ThreadsWithCounter  threads(NTHREADS);
    int                 L               =   cover.getMaxLevelPresent();
    const               Point** ptarr   =   new const Point*[k*NY];
    REAL* distances                     =   new REAL[k*NY];
    FindNearestData     findnearestdata(&threads,*vectorsY,k,L,ptarr,distances);
    
    vectorsY->reset();
//    cout << "\n findNearest";
    cover.findNearest(*vectorsY,findnearestdata,ptarr,distances);
//    cout << "...done findNearest";
    
//    for(int j=0;j<NY;j++) {
//        for(int i=0;i<k;i++){
//            const Point* foundpoint=ptarr[j*k+i];
//            mexPrintf("\n (%d,%d) %u",j,i,foundpoint);
//            if(foundpoint) {
//                mexPrintf("\n %d, %f", vectorsX->getIndex(foundpoint), distances[j*k+i] );
//            } else {
//                mexPrintf("\n -1,-1.0");
//            }
//        }
//    }
    
//    bool test=cover.checkFindNearest(vectorsY,ptarr,k);
//    for(int i=0;i<NY;i++) {
//        if(test==false) {
//            cout << "checkFindNearest failed" << endl;
//        }
//    }
    
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateStructMatrix(1, 1, 2, fnames_out);
    
    mxArray* fout;
    
    dims[0] = k;
    dims[1] = NY;
    fout    = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
    indices = (int*)mxGetData(fout);
    mxSetField(plhs[0],0,fnames_out[0], fout);

    dims[0]=k;
    dims[1]=NY;
#ifdef DOUBLE
    fout=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
#else
    fout=mxCreateNumericArray(ndims,dims,mxSINGLE_CLASS,mxREAL);
#endif
    REAL* outdistances  =   (REAL*)mxGetData(fout);
    mxSetField(plhs[0],0,fnames_out[1], fout);
    
    //mexPrintf("\nDisplaying results...");
    for(int j=0;j<NY;j++) {
        //mexPrintf("\n j=%d", j);
        //const Point* querypoint=ptarr[j*(k+1)];
        //int queryindex=vectorsY.getIndex(querypoint);
        for(int i=0;i<k;i++){
            //mexPrintf("\t i=%d", i);
            const Point* foundpoint=ptarr[j*k+i];
            if(foundpoint) {
                indices[j*k+i]=vectorsX->getIndex(foundpoint);
                outdistances[j*k+i]=distances[j*k+i];
            } else {
                indices[j*k+i]=-1;
                outdistances[j*k+i]=-1.0;
            }
            //mexPrintf("\n %d %d,%f",j,indices[j*k+i],outdistances[j*k+i]);
        }
    }
    
    delete cover;
    delete vectorsY;
    delete vectorsX;
    delete seglist;
    delete [] ptarr;
    delete [] distances;
    
    //mxFree(ptarr);
    //mxFree(distances);    
}

#include "FastSeg.C"
template class SegList<DLPtrListNode<CoverNode> >;










//
//
//// Construct vectors X
//REAL* X=(REAL*)mxGetData(tmp);
//switch (VECTOR_CLASS)   {
//    case VECTORS:
//        *vectorsX = new Vectors(X,NX,dim,DISTANCE_FCN);
//        break;
//    case IMAGES:
//        *vectorsX = new Images(X,NX,sqrt(dim),sqrt(dim),DISTANCE_FCN);
//        break;
//    case MOLECULARSTATES:
//        if( dim % 3==0 )    {
//            *vectorsX = new MolecularStates(X,NX,dim/3);
//        }   else    {
//            mexErrMsgTxt("Number of dimensions should be a multiple of 3\n");
//        }
//        break;
//    default:
//        mexErrMsgTxt("\n Invalid classname.");
//        break;
//}
//
//// The third argument,   withinradius, numfindlevels
//int nfields = mxGetNumberOfFields(prhs[2]);
//if(nfields!=2)                                                                          { mexErrMsgTxt("Third argument should have two fields."); }
//
//tmp=mxGetField(prhs[2],0,within_in[0]);
//if(!mxIsClass(tmp,mxIsClassName))                                                       { mexErrMsgTxt("First field of third argument must be of correct real type\n"); }
//mwSize ndims_radius=mxGetNumberOfDimensions(tmp);
//const mwSize* dims_radius=mxGetDimensions(tmp);
//if(ndims_radius!=2)                                                                     { mexErrMsgTxt("First field of third argument should have 2 dims\n"); }
//if(dims_radius[1]!=1)                                                                   { mexErrMsgTxt("First field of third argument should have dims[1]=1\n"); }
//REAL* pwithinradius=(REAL*)mxGetData(tmp);
//
//tmp=mxGetField(prhs[2],0,within_in[1]);
//if(!mxIsClass(tmp,"int32"))                                                             { mexErrMsgTxt("Third field of fourth argument must be int32\n"); }
//mwSize ndims_findlevels=mxGetNumberOfDimensions(tmp);
//const mwSize* dims_findlevels=mxGetDimensions(tmp);
//if(ndims_findlevels!=2)                                                                 { mexErrMsgTxt("Second field of third argument should have 2 dims\n"); }
//if(dims_findlevels[1]!=1)                                                               { mexErrMsgTxt("Second field of third argument should have dims[1]=1\n"); }
//int* pnumfindlevels=(int*)mxGetData(tmp);
//if(dims_radius[1]!=dims_findlevels[1])                                                  { mexErrMsgTxt("Size mismatch in third argument\n"); }
//
//// Construct vectors Y
//tmp=prhs[3];
//mwSize ndims_Y=mxGetNumberOfDimensions(tmp);
//const mwSize* dims_Y=mxGetDimensions(tmp);
//if(!mxIsClass(tmp,mxIsClassName))                                                       { mexErrMsgTxt("Fourth argument should be of correct real type\n"); }
//int NY=dims_Y[ndims_Y-1];
//int dimY=1;
//for(int i=0;i<ndims_Y-1;i++) {
//    dimY*=dims_Y[i];
//}
//if(dimY!=dim)                                                                           { mexErrMsgTxt("Dimension mismatch\n"); }
//
//REAL* Y=(REAL*)mxGetData(tmp);
//switch (VECTOR_CLASS)   {
//    case VECTORS:
//        *vectorsY = new Vectors(Y,NY,dim,DISTANCE_FCN);
//        break;
//    case IMAGES:
//        *vectorsY = new Images(Y,NY,sqrt(dim),sqrt(dim),DISTANCE_FCN);
//        break;
//    case MOLECULARSTATES:
//        if( dim % 3==0 )    {
//            *vectorsY = new MolecularStates(Y,NY,dim/3);
//        }   else    {
//            mexErrMsgTxt("Number of dimensions should be a multiple of 3\n");
//        }
//        break;
//    default:
//        mexErrMsgTxt("\n Invalid classname.");
//        break;
//}
//
///* Check that third and fourth arguments are compatible */
//if(dims_radius[1]!=1) {
//    if((dims_radius[1]!=NY)||(dims_findlevels[1]!=NY))                                  { mexErrMsgTxt("Mismatch between third and fourth arguments\n"); }
//}
//
//// Get number of threads
//tmp=prhs[4];
//ndims_in=mxGetNumberOfDimensions(tmp);
//dims_in=mxGetDimensions(tmp);
//if(!mxIsClass(tmp,"int32"))                                                             { mexErrMsgTxt("Fifth argument should be int32\n"); }
//*NTHREADS=*(int*)mxGetData(tmp);
//
//// Get number of trees
//tmp=prhs[5];
//ndims_in=mxGetNumberOfDimensions(tmp);
//dims_in=mxGetDimensions(tmp);
//if(!mxIsClass(tmp,"int32"))                                                             { mexErrMsgTxt("Sixth argument should be int32\n"); }
//*NTREES=*(int*)mxGetData(tmp);
//

