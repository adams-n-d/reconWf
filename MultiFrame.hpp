//#include "script.hpp"

int GetCompositeImageSize(int *mLarge, int *nLarge, const int m, const int n, const double deltax, const double deltay,const int nframes)
{

    double *Tx = new double[nframes];
    double *Ty = new double[nframes];

    int err;

    err = GenerateMotionData(Tx,Ty,nframes,deltax,deltay);
    if(err) cout << "Failure to generate the motion data Tx and Ty" << endl;

    // generate coordinates
    double *X = new double[m*n*nframes];
    double *Y = new double[m*n*nframes];
    err = GetPixelCenter2D(X,Y,m,n);
    if(err) cout << "Failure to generate the coordinate system" << endl;

    // shift the coordinates
    for(int i=1; i<nframes; ++i) {
        err = TransformCoordinates2D(&X[i*m*n], &Y[i*m*n], Tx[i], Ty[i], X, Y, m, n);
        if(err) cout << "Failure to shift " << i+1 << "frame coordinates" << endl;
    }

    // normalize the shifted coordinates
    double minX = n, minY = n, maxX = 0, maxY = 0;
    for(int i=0; i<m*n*nframes; ++i ) {
        minX = (minX < X[i]) ? minX : X[i];
        minY = (minY < Y[i]) ? minY : Y[i];
        maxX = (maxX < X[i]) ? X[i] : maxX;
        maxY = (maxY < Y[i]) ? Y[i] : maxY;
    }

    //int mLarge, nLarge;
    *mLarge = static_cast<int>(round(maxY-minY+1) > round(maxX-minX+1) ? round(maxY-minY+1) : round(maxX-minX+1));

    if(*mLarge%2 != 0) *mLarge += 1;
    *nLarge = *mLarge;

    delete Tx;
    delete Ty;
    delete X;
    delete Y;

    return(EXIT_SUCCESS);
}

int ConstructMotionMatrixIndices(int *rowIdx, int *colIdx, double *vals, const double *vecI, const double *vecJ, const int m, const int n, const int nframes, const int *rowIdxLength)

{        

    // number of nonzeros per row
    int numEntries = 4;
    int err;

    // find the starting point of rowIdx, colIdx and vals for each frame
    int *start = new int[nframes];
    start[0] = 0;
    int preRowLength = rowIdxLength[0];
    for (int i=1; i<nframes; ++i) {
        start[i] = preRowLength;
        preRowLength += rowIdxLength[i];
    }

    for (int i=0; i<nframes; ++i) {
       err = BuildInterpMatrix(&vecI[i*m*n], &vecJ[i*m*n], &rowIdx[start[i]], &colIdx[start[i]*numEntries], &vals[start[i]*numEntries], m, n, rowIdxLength[i], i);
       for(int j=0; j<rowIdxLength[i]; ++j) {
           rowIdx[start[i]+j] += i*m*n;
       }
    }

   delete start;

   return(EXIT_SUCCESS);
}



int ConstructSparseDiagonalMatrixIndices(int *iwLarge, int *jwLarge, double *valsLarge, const int *iw, const int *jw, const double *vals, int nnz, int nframes, int m, int n)
{
    // m is the number of rows in each small sparse R or Wp.
    // n is the number of columns in each small sparse R or Wp.
    // nnz is the number of non-zeros per row in each small sparse R or Wp.
    // int *iwLarge = new int[nframes*nnz];
    // int *jwLarge = new int[nframes*nnz];
    // double *valsLarge = new double[nframes*nnz];
    
    for(int i=0; i<nframes; ++i)
        for(int j=0; j<nnz; ++j) {
            iwLarge[i*nnz+j] = iw[j] + i*m;
            jwLarge[i*nnz+j] = jw[j] + i*n;
            valsLarge[i*nnz+j] = vals[j];
        }

    return(EXIT_SUCCESS);

}


