//Compile command: mex -O -DDOUBLE findwithin.C ThreadsWithCounter.C IDLList.C IDLListNode.C Vector.C Cover.C findWithin.C Point.C CoverNode.C EnlargeData.C FindWithinData.C -lpthread

/* findwithin take five arguments:
 The first argument is what was returned by covertree
 The second argument is the array that was passed as the second argument to
 covertree
 The third argument is a struct whose first field is "distances"
 of class double
 and whose second field is "numlevels" of class int32.
 both are either 1x1 or NYx1 where NY is computed from the fourth argument
 as below.
 The fourth argument is the array of search vectors.
 The fifth argument, of class int32, is the number of threads (0 for serial).
 
 findwithin returns a struct whose first member is "pi" of class int32
 and is NYx2; whose second members "indices" is of class int32 and is
 totalfound x 1; and whose third member "distances" is of class double
 and is totalfound x 1. The (i,1) entry of pi is the number found
 correspoinding the j'th search vector and the (i,2) entry of indices is the
 offset in found of the indices of the corresponding found vectors with
 respect to the second argument. The distances corresponding to a given
 query point are sorted in nondecreasing order.
 */

#include <math.h>
#include "mex.h"
#include "Cover.H"
#include "CoverForest.H"
#include "Vector.H"
#include "ThreadsWithCounter.H"
#include "FindWithinData.H"
#include "Distances.H"
#include "VectorClassNames.H"


const char* fnames_in[]={
    "theta",
    "outparams",
    "radii",
    "levels",
    "ncallstogetdist",
};

const char* within_in[]={
    "distances",
    "numlevels"
};

const char* fnames_out[]={
    "pi",
    "indices",
    "distances"
};


int     dim=0;
REAL    theta=0.0;
REAL    mu=0.0;
mwSize  dims[2];
int     ndims=2;


void ConvertMatStructToCoverTree( const mxArray *prhs[], Distance_Mode* DISTANCE_FCN, VectorClassNames* VECTOR_CLASS, Vectors** vectorsX, Vectors** vectorsY, int* minlevel, int* numlevels, int*q,
                                 unsigned int *NTHREADS, unsigned int *NTREES );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    Distance_Mode DISTANCE_FCN;
    VectorClassNames VECTOR_CLASS;
    
#ifdef DOUBLE
    char mxIsClassName[] = "double";
#else
    char mxIsClassName [] = "single";
#endif
    
    // Check for proper number of arguments.
    if (nrhs != 6)                                                                          { mexErrMsgTxt("Six inputs required."); }
    // Perform range search with cover trees or cover forests
    if( NTREES<=1 )  {
        Distance_Mode DISTANCE_FCN;
        VectorClassNames VECTOR_CLASS;
        Vectors* vectorsX, vectorsY;
        int minlevel, numlevels, q;
        unsigned int NTHREADS;
        unsigned int NTREES;
        
        // Re-create parameters of hte cover tree
        ConvertMatStructToCoverTree( prhs, &DISTANCE_FCN, &VECTOR_CLASS, &vectorsX, &vectorsY, &minlevel, &numlevels, &q, &NTHREADS, &NTREES );
        
        // Re-create covertree
        CoverIndices coverindices(theta,numlevels,minlevel,NX,q);
        SegList<DLPtrListNode<CoverNode> > seglist(NX);
        Cover cover(*vectorsX,seglist,coverindices,DISTANCE_FCN);
        
        Cover::DescendList* descendlists=new Cover::DescendList[NY];
        int totalfound=0;
        ThreadsWithCounter threads(NTHREADS);
        if(dims_radius[0]==1) {
            FindWithinData findwithindata(&threads,*vectorsY,*pwithinradius,*pnumfindlevels,descendlists);
            totalfound=cover.findWithin(*vectorsY,findwithindata,descendlists);
        } else {
            FindWithinData findwithindata(&threads,*vectorsY,pwithinradius,pnumfindlevels,descendlists);
            totalfound=cover.findWithin(*vectorsY,findwithindata,descendlists);
        }
        vectorsY->reset();
        
        // Create matrix for the return argument.
        plhs[0] = mxCreateStructMatrix(1, 1, 3, fnames_out);
        
        mxArray* fout;
        
        dims[0] = NY;
        dims[1] = 2;
        fout    = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
        int* pi = (int*)mxGetData(fout);  //numfound and offsets
        
        mxSetField(plhs[0],0,fnames_out[0],fout);
        
        dims[0] = totalfound;
        dims[1] = 1;
        fout    = mxCreateNumericArray(ndims,dims,mxUINT32_CLASS,mxREAL);
        INDEX* indices=(INDEX*)mxGetData(fout);
        
        mxSetField(plhs[0],0,fnames_out[1], fout);
        
        dims[0]=totalfound;
        dims[1]=1;
#ifdef DOUBLE
        fout =mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
#else
        fout =mxCreateNumericArray(ndims,dims,mxSINGLE_CLASS,mxREAL);
#endif
        REAL* distances=(REAL*)mxGetData(fout);
        
        mxSetField(plhs[0],0,fnames_out[2], fout);
        
        // Fill in outputs with flattened lists of neighbors
        FindWithinflattenDescendLists((INDEX)totalfound,(INDEX *)indices,distances,(INDEX)NY,(INDEX *)pi,*vectorsX,descendlists);
        //bool test=checkFlattenDescendList(totalfound,indices,NY,pi,vectorsX,descendlists);
        
        delete [] descendlists;
        
    } else {
        CoverIndices coverindices(theta,numlevels,minlevel,NX,q);
        SegList<DLPtrListNode<CoverNode> > seglist(NX);
        CoverForest coverForest(*vectorsX,seglist,coverindices,DISTANCE_FCN);
        
        int totalfound=0;
        if(dims_radius[0]==1) {
            INDEX *pi_forest, *idxarr_forest;
            REAL *distances_forest;
            
            totalfound = coverForest.findWithin( *vectorsY,*pwithinradius,&pi_forest, &idxarr_forest,&distances_forest );
            
            // Create matrix for the return argument.
            plhs[0] = mxCreateStructMatrix(1, 1, 3, fnames_out);
            
            mxArray* fout;
            
            dims[0] = NY;
            dims[1] = 2;
            fout    = mxCreateNumericArray(ndims,dims,mxINT32_CLASS,mxREAL);
            INDEX* pi = (INDEX*)mxGetData(fout);  //numfound and offsets
            memcpy( pi, pi_forest, 2*NY*sizeof(INDEX) );
            
            mxSetField(plhs[0],0,fnames_out[0],fout);
            
            dims[0] = totalfound;
            dims[1] = 1;
            fout    = mxCreateNumericArray(ndims,dims,mxUINT32_CLASS,mxREAL);
            INDEX* indices=(INDEX*)mxGetData(fout);
            memcpy( indices, idxarr_forest, totalfound*sizeof(INDEX) );
            
            mxSetField(plhs[0],0,fnames_out[1], fout);
            
            dims[0]=totalfound;
            dims[1]=1;
#ifdef DOUBLE
            fout =mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
#else
            fout =mxCreateNumericArray(ndims,dims,mxSINGLE_CLASS,mxREAL);
#endif
            REAL* distances=(REAL*)mxGetData(fout);
            memcpy( distances, distances_forest, totalfound*sizeof(REAL) );
            
            mxSetField(plhs[0],0,fnames_out[2], fout);
            
            delete [] pi_forest;
            delete [] idxarr_forest;
            delete [] distances_forest;
            
        } else {
            mexErrMsgTxt("Searches with multiple radii not implemented with multiple trees.\n");
        }
    }
    
    delete vectorsY;
    delete vectorsX;
    
}










#include "FastSeg.C"
template class SegList<DLPtrListNode<CoverNode> >;
