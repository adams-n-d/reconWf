#define NOCOL 128 
#define NFRAMES 1
#include <iostream>
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "EpetraExt_CrsMatrixIn.h"

#include "AztecOO.h"

//#include <Epetra_Import.h>
//#include <Epetra_LocalMap.h>
//#include <Epetra_CrsGraph.h>
//#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_BlockMapIn.h>
#include <EpetraExt_CrsMatrixIn.h>
//#include <EpetraExt_VectorIn.h>
#include "Epetra_MultiVector.h"
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_MultiVectorOut.h>

//#include "BelosLinearProblem.hpp"
//#include "BelosEpetraAdapter.hpp"
//#include "BelosLSQRSolMgr.hpp"

using namespace EpetraExt;

int Vectorize(const double image[][NOCOL], double *vec, const int m, const int n)
{
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            vec[i+j*n] = image[i][j];
 
    return(EXIT_SUCCESS);
}

int Reshape(const double *vec, double image[][NOCOL], const int m, const int n)
{
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            image[i][j] = vec[i+j*n];

    return(EXIT_SUCCESS);
}

int GetPixelCenter2D(double *X, double *Y, const int m, const int n)
{
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j) {
            X[i+j*n] = j*1.0;
            Y[i+j*n] = (m-i-1)*1.0;
        }

    return(EXIT_SUCCESS);
}

int GenerateMotionData(double *Tx, double *Ty, const int nframes, const double deltax, const double deltay)
{
    for(int i=0; i<nframes; ++i) {
        Tx[i] = i*deltax;
        Ty[i] = i*deltay;
    }

    return(EXIT_SUCCESS);
}

int TransformCoordinates2D(double *Xnew, double *Ynew, const double Tx, const double Ty, const double *X, const double *Y, const int m, const int n)
{

    for(int i=0; i<m*n; ++i) {
        Xnew[i] = X[i] + Tx;
        Ynew[i] = Y[i] -Ty;
    }

    return(EXIT_SUCCESS);
}

int SpaceToCoordinates2D(const double *X, const double *Y, double *I, double *J, const int m, const int n)
{
    for(int i=0; i<m*n; ++i) {
        I[i] = m-1-Y[i];
        J[i] = X[i];
    }

    return(EXIT_SUCCESS);
}

int BuildInterpMatrix(const double *vecI, const double *vecJ, Epetra_Map map, Epetra_CrsMatrix *A, const int m, const int n)
{
    int myDimensions = map.NumMyElements();
    int *myGlobalDimensions = map.MyGlobalElements();
    int numEntries = 4;

    int *i0 = new int[m*n];
    int *i1 = new int[m*n];
    int *j0 = new int[m*n];
    int *j1 = new int[m*n];

    // 
    for(int index = 0; index < m*n; ++index) {
        i0[index] = static_cast<int>(floor(vecI[index]));
        i1[index] = i0[index] +1;
        j0[index] = static_cast<int>(floor(vecJ[index]));
        j1[index] = j0[index] +1;
    }

    //
    int rowIdxLength = 0;
    for (int index = 0; index < m*n; ++index) {
        if (i0[index]>=0 && i1[index]<m && j0[index] >= 0 && j1[index] < n) {
            rowIdxLength += 1;
        }
    }

    //
    int * rowIdx = new int [rowIdxLength];
    int indexForRowIdx = 0;
    
    for (int index = 0; index < m*n; ++index) {
        if (i0[index]>=0 && i1[index]<m && j0[index] >= 0 && j1[index] < n) {
            rowIdx[indexForRowIdx] = index;
            indexForRowIdx += 1;
        }
    }

    // compute the corresponding column index
    int *columnIdx1 = new int[rowIdxLength];
    int *columnIdx2 = new int[rowIdxLength];  
    int *columnIdx3 = new int[rowIdxLength];
    int *columnIdx4 = new int[rowIdxLength];
    
    for (int index = 0; index < rowIdxLength; ++index) {
        columnIdx1[index] = i0[rowIdx[index]] + m*(j0[rowIdx[index]]);
        columnIdx2[index] = i1[rowIdx[index]] + m*(j0[rowIdx[index]]);
        columnIdx3[index] = i0[rowIdx[index]] + m*(j1[rowIdx[index]]);
        columnIdx4[index] = i1[rowIdx[index]] + m*(j1[rowIdx[index]]);
    }

    // compute the weighting values
    double deltai, deltaj;
    double * weight1, *weight2, *weight3, *weight4;
    weight1 = new double [rowIdxLength];
    weight2 = new double [rowIdxLength];
    weight3 = new double [rowIdxLength];
    weight4 = new double [rowIdxLength];

    for (int index = 0; index < rowIdxLength; ++index) {
        deltai = vecI[rowIdx[index]] - i0[rowIdx[index]];
        deltaj = vecJ[rowIdx[index]] - j0[rowIdx[index]];
        weight1[index] = (1-deltai) * (1-deltaj);
        weight2[index] = (1-deltaj) * deltai;
        weight3[index] = (1-deltai) * deltaj;
        weight4[index] = (deltai) * deltaj;
    }

    // fill in the matrix data
    double *values = new double[4];
    int *indices = new int[4];
    for (int index = 0; index < myDimensions; ++index) {
        for (int indexForRowIdx = 0; indexForRowIdx < rowIdxLength; ++indexForRowIdx) {
            if (myGlobalDimensions[index] == rowIdx[indexForRowIdx]) {
                values[0] = weight1[indexForRowIdx];
                values[1] = weight2[indexForRowIdx];
                values[2] = weight3[indexForRowIdx];
                values[3] = weight4[indexForRowIdx];
                
                // global column indices
                indices[0] = columnIdx1[indexForRowIdx];
                indices[1] = columnIdx2[indexForRowIdx];
                indices[2] = columnIdx3[indexForRowIdx];
                indices[3] = columnIdx4[indexForRowIdx];
                (*A).InsertGlobalValues(myGlobalDimensions[index], numEntries, values, indices);
            }
        }
    }
    (*A).FillComplete();

    //  free memory;
    //delete vecI;
    //delete vecJ;
    delete i0;
    delete i1;
    delete j0;
    delete j1;

    delete rowIdx;
    delete weight1;
    delete weight2;
    delete weight3;
    delete weight4;

    delete columnIdx1;
    delete columnIdx2;
    delete columnIdx3;
    delete columnIdx4;
  
    delete values;
    delete indices;

    return(EXIT_SUCCESS);
}


int BuildInterpMatrix(const double *vecI, const double *vecJ, int *rowIdx, int *columnIdx, double *weight, const int m, const int n, const int rowIdxLength, const int frameIdx)
{
    
    int numEntries = 4;

    int *i0 = new int[m*n];
    int *i1 = new int[m*n];
    int *j0 = new int[m*n];
    int *j1 = new int[m*n];

    // 
    for(int index = 0; index < m*n; ++index) {
        i0[index] = static_cast<int>(floor(vecI[index]));
        i1[index] = i0[index] +1;
        j0[index] = static_cast<int>(floor(vecJ[index]));
        j1[index] = j0[index] +1;
    }


    //
    int indexForRowIdx = 0;
    for (int index = 0; index < m*n; ++index) {
        if (i0[index]>=0 && i1[index]<m && j0[index] >= 0 && j1[index] < n) {
            rowIdx[indexForRowIdx] = index;
            indexForRowIdx += 1;
        }
    }

    // compute the corresponding column index
    //int *columnIdx = new int[rowIdxLength*numEntries];
    
    for (int index = 0; index < rowIdxLength; ++index) {
        columnIdx[index*numEntries] = i0[rowIdx[index]] + m*(j0[rowIdx[index]]);
        columnIdx[index*numEntries +1] = i1[rowIdx[index]] + m*(j0[rowIdx[index]]);
        columnIdx[index*numEntries +2] = i0[rowIdx[index]] + m*(j1[rowIdx[index]]);
        columnIdx[index*numEntries +3] = i1[rowIdx[index]] + m*(j1[rowIdx[index]]);
    }

    // compute the weighting values
    double deltai, deltaj;
    //double * weight;
    //weight = new double [rowIdxLength*numEntries];

    for (int index = 0; index < rowIdxLength; ++index) {
        deltai = vecI[rowIdx[index]] - i0[rowIdx[index]];
        deltaj = vecJ[rowIdx[index]] - j0[rowIdx[index]];
        weight[index*numEntries] = (1-deltai) * (1-deltaj);
        weight[index*numEntries +1] = (1-deltaj) * deltai;
        weight[index*numEntries +2] = (1-deltai) * deltaj;
        weight[index*numEntries +3] = (deltai) * deltaj;
    }

    //  free memory;
    delete i0;
    delete i1;
    delete j0;
    delete j1;

    return(EXIT_SUCCESS);
}

int FindInterpMatrixLength(const double *vecI, const double *vecJ, const int m, const int n, int *rowIdxLength)
{
    
    int *i0 = new int[m*n];
    int *i1 = new int[m*n];
    int *j0 = new int[m*n];
    int *j1 = new int[m*n];

    // 
    for(int index = 0; index < m*n; ++index) {
        i0[index] = static_cast<int>(floor(vecI[index]));
        i1[index] = i0[index] +1;
        j0[index] = static_cast<int>(floor(vecJ[index]));
        j1[index] = j0[index] +1;
    }

    //
    *rowIdxLength = 0;
    for (int index = 0; index < m*n; ++index) {
        if (i0[index]>=0 && i1[index]<m && j0[index] >= 0 && j1[index] < n) {
            *rowIdxLength += 1;
        }
    }

    delete i0;
    delete i1;
    delete j0;
    delete j1;

    return(EXIT_SUCCESS);
}

int rest2D(int *iw, int *jw, double *vals, int *nnz, const int m, const int n, const int ssp)
{
    if(n%ssp != 0 || m%ssp != 0) {
        cout << "m or n must be multiple of ssp!" << endl;
        return(1);
    }
    
    *nnz = m*n;
    int idx = 0;

    for(int l=0; l<n/ssp; ++l)
        for(int i=0; i<m/ssp; ++i)
            for(int j=0; j<ssp; ++j)
                for(int k=0; k<ssp; ++k) {
                    jw[idx] = l*m*ssp+i*ssp+j*m+k;
                    iw[idx] = l*m/ssp+i;
                    vals[idx] = 1.0/(ssp*ssp);
                    idx += 1;
                }


    return(EXIT_SUCCESS);
}

int MakeMask(const int nf, const double r1, double *mask)
{
    double h = 2.0/nf;
    double *X1 = new double[nf*nf];
    double *X2 = new double[nf*nf];

    for(int i=0; i<nf; ++i)
        for(int j=0; j<nf; ++j) {
            X1[i+nf*j] = -1+i*h;
            X2[i+nf*j] = -1+j*h;
        }

    for(int i=0; i<nf*nf; ++i)
        mask[i] = sqrt(X1[i]*X1[i] + X2[i]*X2[i]);
    
    for(int i=0; i<nf*nf; ++i) {
        if(mask[i] <= r1) mask[i] = 1.0;
        else mask[i] = 0.0;
    }

    delete X1;
    delete X2;

    return(EXIT_SUCCESS);
}

int ConstructPupilWindowMatrixIndices(int *iw, int *jw, double *mask, const int mf, const int nf, const int mLarge, const int nLarge,const double r1)
{
    int err;
    err = MakeMask(nf,r1,mask);

    int *jw1 = new int[mf*nf];
    int *jw2 = new int[mf*nf];
    
    for(int i=0; i<nf; ++i)
    for(int j=0; j<mf; ++j) {
        jw1[i*mf+j] = mLarge*i;
        jw2[i*mf+j] = j;
    }

    for(int i=0; i<mf*nf; ++i) {
        iw[i] = i;
        jw[i] = jw1[i] + jw2[i];
    }

    delete jw1;
    delete jw2;

    return(EXIT_SUCCESS);
}


