#include <matio.h>
#include <Epetra_MultiVector.h>
using namespace EpetraExt;


// Read matrices and vectors from a .mat file with multiple variables
int readMat(const char *fileList, char* varName, Epetra_MultiVector *& nGradMeas, Epetra_Map &map){

    mat_t *mat;
    mat = Mat_Open(fileList, MAT_ACC_RDONLY);

    matvar_t *ngrad = Mat_VarReadInfo(mat, varName);
    Mat_VarPrint(ngrad, 1);

    int m, n, z;
    m = ngrad->dims[0];
    n = ngrad->dims[1];    
    //z = ngrad->dims[2];
    z = 1;
    int start[3]={0,0,0},stride[3]={1,1,1},edge[3]={m,n,z};
    void *ndat = malloc(m*n*z*sizeof(double));
    int err = Mat_VarReadData(mat, ngrad, ndat, start, stride, edge);
    double* ndat_dbl = (double*)ndat;
    int i,j;
    /*for(i=0; i < m*n; i++){
      cout << "t " << ndat_dbl[i] << endl;
    }*/
    //Epetra_MultiVector gm(Copy, map, xdat_dbl, m, n);
    //xGradMeas = &gm;
    //cout << gm << endl;   

    /*double** epetArr = (double**)malloc(n*sizeof(double*));
    for(i=0; i<n; i++){
        epetArr[i] = &ndat_dbl[i*n];
    }*/
//    xGradMeas = new Epetra_MultiVector(Copy, map, epetArr, n); 
    int one = 1;
    Epetra_MultiVector GradMeas(Copy, map, ndat_dbl, m*n*z, one);
   // xGradMeas = new Epetra_MultiVector(map, n);
    nGradMeas = new Epetra_MultiVector(GradMeas);
//    Epetra_MultiVector xGradMeas(*GradMeas);
/*    double ** test;
    (*xGradMeas).ExtractView(&test);
    for(i=0; i < m; i++){
      for(j=0; j < n; j++){
        cout << test[i][j] << endl;
      }
    }
  */  //cout << (*xGradMeas) << endl;

 //   xGradMeas->Epetra_MultiVector(Copy, map, epetArr, n); 

    return(EXIT_SUCCESS);
}

// Read matrices and vectors from a .mat file with multiple variables
int readMat(const char *fileList, char* varName, Epetra_CrsMatrix *& spMatrix, Epetra_Map &map){

    mat_t *mat;
    mat = Mat_Open(fileList, MAT_ACC_RDONLY);

    matvar_t *ngrad = Mat_VarReadInfo(mat, varName);
    Mat_VarPrint(ngrad, 1);

    int m, n, z;
    m = ngrad->dims[0];
    n = ngrad->dims[1];    
    //z = ngrad->dims[2];
    z = 1;
    int start[3]={0,0,0},stride[3]={1,1,1},edge[3]={m,n,z};
    void *ndat = malloc(m*n*z*sizeof(double));
    int err = Mat_VarReadData(mat, ngrad, ndat, start, stride, edge);
    double* ndat_dbl = (double*)ndat;
    int i,j;
    /*for(i=0; i < m*n; i++){
      cout << "t" << xdat_dbl[i] << endl;
    }*/
    //Epetra_MultiVector gm(Copy, map, xdat_dbl, m, n);
    //xGradMeas = &gm;
    //cout << gm << endl;   

    /*double** epetArr = (double**)malloc(n*sizeof(double*));
    for(i=0; i<n; i++){
        epetArr[i] = &ndat_dbl[i*n];
    }*/
//    xGradMeas = new Epetra_MultiVector(Copy, map, epetArr, n); 
    int one = 1;
    //Epetra_CrsMatrix spMat(Copy, map, ndat_dbl, m*n*z, one);
   // xGradMeas = new Epetra_MultiVector(map, n);
    double *vals = new double[n];
    int *indc = new int[n];
    for(i=0; i < n; ++i){
       indc[i]=i;
    }


    spMatrix = new Epetra_CrsMatrix(Copy, map, n/2);
    for(i=0; i < m; ++i){
       for(j=0; j < n; ++j){
cout << "v: " << ndat_dbl[n*i +j] << endl;
          vals[j] = ndat_dbl[n*i + j];
cout << "pupils: " << i << ", " << j << ": " <<vals[j]<<endl;;
       }
       (*spMatrix).InsertGlobalValues(i, n, vals, indc);

    }
    (*spMatrix).FillComplete(map, map);
//    Epetra_MultiVector xGradMeas(*GradMeas);
/*    double ** test;
    (*xGradMeas).ExtractView(&test);
    for(i=0; i < m; i++){
      for(j=0; j < n; j++){
        cout << test[i][j] << endl;
      }
    }
  */  //cout << (*xGradMeas) << endl;

 //   xGradMeas->Epetra_MultiVector(Copy, map, epetArr, n); 

    return(EXIT_SUCCESS);
}

// Read matrices and vectors.
int readMapMatrix(const char *fileList[], const int numFile, Epetra_Comm &comm,
		  Epetra_MultiVector *& xGradMeas, Epetra_MultiVector *& yGradMeas,Epetra_Map &map) 
{
    int err;

    // read in the measured x- and y- gradients (on coarse grids)
    err = MatrixMarketFileToMultiVector(fileList[0], map, xGradMeas);
    if (err != 0) {
        cout << "Failure reading the x- gradients" << endl;
        return err;
    }

    err = MatrixMarketFileToMultiVector(fileList[1], map, yGradMeas);
    if (err != 0) {
        cout << "Railure reading the y- gradients" << endl;
        return err;
    }

    return(EXIT_SUCCESS);
}

// This is a Trilinos function, but not accessible
// from EpetraExt_BlockMapIn.h.
int MatrixMarketFileToRowMap(const char* filename,
                             const Epetra_Comm& comm,
                              Epetra_BlockMap*& rowmap)
{
    FILE* infile = fopen(filename, "r");
    MM_typecode matcode;
 
    int err = mm_read_banner(infile, &matcode);
    if (err != 0) return(err);
 
    if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) ||
        !mm_is_real(matcode)   || !mm_is_general(matcode)) {
       return(-1);
    }
 
    int numrows, numcols;
    err = mm_read_mtx_array_size(infile, &numrows, &numcols);
    if (err != 0) return(err);
 
    fclose(infile);
 
    rowmap = new Epetra_BlockMap(numrows, 1, 0, comm);
    return(0);
 }


