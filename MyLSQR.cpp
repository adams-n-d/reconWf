#include "MatVec.h"
#include "script.hpp"
#include "MultiFrame.hpp"
#include "ReadIn.h"
#include <EpetraExt_MatrixMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Time.h>

using namespace EpetraExt;
bool verbose;
bool lead_proc;
int MyPID;


//Construct the Fried/Hudgens geometry matrix
int ConstructWFOperator(int n, Epetra_FECrsMatrix *& A, Epetra_MpiComm &comm, Epetra_CrsMatrix *& Wp){
    int i, j, k;
    int err, err2;

    int n2 = n*n;
    int idx[4];

    Epetra_Map rowMap(2*n2, 0, comm);
    Epetra_Map colMap(n2, 0, comm);

    Epetra_Map DrowMap(n2, 0, comm);
    Epetra_Map DcolMap(n2, 0, comm);

    Epetra_FECrsMatrix *P = new Epetra_FECrsMatrix(Copy, DrowMap, 1);
    Epetra_CrsMatrix *Px = new Epetra_FECrsMatrix(Copy, rowMap, 1);
    Epetra_CrsMatrix *Py = new Epetra_FECrsMatrix(Copy, rowMap, 1);

    //put n*n pupil mask into n2*n2 diagonal matrix
    int cnt =0;
    k=0;
    int one = 1;
    double pval[1];
    int pidx[1];

    int nZWp;
    double *pupVals = new double[n];
    int *pupIdx = new int[n];


    for(i=0; i < n; i++){
    if( (*Wp).MyGRID(i) ){  
      (*Wp).ExtractGlobalRowCopy(i, n, nZWp, pupVals, pupIdx);
      for(j=0; j < nZWp; j++){
        pval[0] = pupVals[j];
	
	k = i + pupIdx[j]*n;
        pidx[0] = k;

        err = (*P).InsertGlobalValues(k, one, pval, pidx);

        if(err !=0){
          cout << "P_ insert error!" << endl;
          cout << "Px error code: " << err << endl;
          //cout << "Py error code: " << err2 << endl;

          MPI_Finalize();
          exit(1);
        }
      }
    }
    }
    (*P).GlobalAssemble();


// Create [P*Dx] and [P*Dy], then stack to get operator
    Epetra_FECrsMatrix *Dx = new Epetra_FECrsMatrix(Copy, DrowMap, 4);
    Epetra_FECrsMatrix *Dy = new Epetra_FECrsMatrix(Copy, DrowMap, 4);
    Epetra_FECrsMatrix *PDx = new Epetra_FECrsMatrix(Copy, DrowMap, 4);
    Epetra_FECrsMatrix *PDy = new Epetra_FECrsMatrix(Copy, DrowMap, 4);


//create Dx and Dy, which represent kron(F, H) and kron(H, F) respectively
// where F and H are the Fried and Hudgens geometry matrices
    double vals1[4] = {-0.5, 0.5, -0.5, 0.5};
    double vals2[4] = {-0.5, -0.5, 0.5, 0.5};

    int numEntries;
    double stubVals1[2] = {2*vals1[0], 2*vals1[1]};
    double stubVals2[2] = {2*vals2[0], 2*vals2[2]};
    int lim = n2 - n;
    for(i=0; i < lim; i++){
        if( (*Dx).MyGRID(i) ){  
            if( (i+1) % n != 0){
	        numEntries = 4;
                idx[0] = i;
                idx[1] = i+1;
                idx[2] = i+n;
                idx[3] = i+n+1;
                (*Dx).InsertGlobalValues(i, numEntries, vals1, idx);
                (*Dy).InsertGlobalValues(i, numEntries, vals2, idx);

	        } else{
                    numEntries = 2;
                    idx[0] = i;
                    idx[1] = i + n;
                    (*Dy).InsertGlobalValues(i, numEntries, stubVals2, idx);
                }

        }
    }
    for(i=lim; i < n2-1; i++){
    if( (*Dx).MyGRID(i) ){  
        numEntries = 2;
        idx[0] = i;
        idx[1] = i + 1;
        (*Dx).InsertGlobalValues(i, numEntries, stubVals1, idx);
    }
    }
    (*Dx).GlobalAssemble();
    (*Dy).GlobalAssemble();

/*     RowMatrixToMatrixMarketFile(  "Data/Dx.mm", *Dx );
     RowMatrixToMatrixMarketFile(  "Data/Dy.mm", *Dy );
*/
    EpetraExt::MatrixMatrix::Multiply(*P, 0, *Dx, 0, *PDx);
    EpetraExt::MatrixMatrix::Multiply(*P, 0, *Dy, 0, *PDy);

    (*PDx).FillComplete();
    (*PDy).FillComplete();
            //RowMatrixToMatrixMarketFile(  "Data/PDx.mm", *(PDx) );
            //RowMatrixToMatrixMarketFile(  "Data/PDy.mm", *(PDy) );

    A = new Epetra_FECrsMatrix(Copy, rowMap, 4);
    int nZPx;
    int nZPy;
    double *xVals = new double[4];
    double *yVals = new double[4];
    int *xIdx = new int[4];
    int *yIdx = new int[4];
    for(i=0; i < n2; i++){
        if( (*Dx).MyGRID(i) ){  

        (*PDx).ExtractGlobalRowCopy(i, n2, nZPx, xVals, xIdx);
        (*PDy).ExtractGlobalRowCopy(i, n2, nZPy, yVals, yIdx);

        (*A).InsertGlobalValues(i, nZPx, xVals, xIdx);
        (*A).InsertGlobalValues(i + n2, nZPy, yVals, yIdx);
        }
    }
    (*A).GlobalAssemble(colMap, rowMap);

    RowMatrixToMatrixMarketFile(  "Data/A.mm", *A );
//this code verifies the assignment and extraction methods are working

//    cout << "crsmatrix A is " << (*A).NumGlobalRows() << " by " << (*A).NumGlobalCols() << ", norm is " << (*A).NumGlobalNonzeros()  << endl;

    return 0;
}


int main(int argc, char **argv)
{


    #ifdef EPETRA_MPI
        MPI_Init(&argc,&argv);
        Epetra_MpiComm comm(MPI_COMM_WORLD);

        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    #else
        Epetra_SerialComm comm;
        int world_size = 1;
    #endif

    //each frame will have frm_sz processors dedicated
    int frm_sz = world_size / NFRAMES;
    // Get our processor ID (0 if running serial)
    MyPID = comm.MyPID();
    // Verbose flag: only processor 0 should print
    verbose = (MyPID==0);
    // The lead_proc is the lowest numbered process for each frame
    lead_proc = (MyPID%frm_sz==0);

    int MyFrame = MyPID / frm_sz;
    int MyLdr = MyFrame*frm_sz; //process in group with lowest rank#

    if(verbose){
        cout << "There are " << NFRAMES << " frames" << endl;
        cout << "Each frame has " << frm_sz << "procs" << endl;
    }
   
    // define the dimensions
    int m = NOCOL;
    int n = NOCOL;
    int nframes = NFRAMES;
    int i, j;

    //create groups by the frame they will work on
    MPI_Comm MPI_frame_comm[NFRAMES]; //a communicator for each frame
    MPI_Group frm_grp[NFRAMES];
    int * fRanks = (int*) malloc(sizeof(int)*frm_sz);
    for(i=0; i < frm_sz; i++){
      fRanks[i] = MyLdr + i;
    }

    MPI_Group glbgrp;
    MPI_Comm_group(MPI_COMM_WORLD, &glbgrp);
    MPI_Group_incl(glbgrp, frm_sz, fRanks, &frm_grp[MyFrame]);
    MPI_Comm_create(MPI_COMM_WORLD, frm_grp[MyFrame], &MPI_frame_comm[MyFrame]);

    //set up communicators
    Epetra_MpiComm frm_comm(MPI_frame_comm[MyFrame]);
    cout << "Process: " << MyPID << "\tFrame: " << MyFrame<< "\tMy comm is " << frm_comm << endl;

    const int * sz = &m;
    int err; // for implementing functions

    Epetra_FECrsMatrix *A = NULL;
    // Declaration of vectors.
    Epetra_MultiVector *xGradMeas = NULL;
    Epetra_MultiVector *yGradMeas = NULL;

    Epetra_CrsMatrix *Wp = NULL;
    const char * pmaskFl = "pm128.txt";
    const char * Wpfile = "Data/Wp.mm";
    char * Wpvar = "p_mask";
    //err = readMat(Wpfile, Wpvar, Wp, Wp_Map);
//cout << "lead proc, my PID: " << MyPID << endl;
    err = MatlabFileToCrsMatrix(pmaskFl, frm_comm, Wp);
    Epetra_Map Wp_Map(n, 0, frm_comm);
//    err = MatrixMarketFileToCrsMatrix(Wpfile, Wp_Map, Wp);
    (*Wp).FillComplete();
    err = RowMatrixToMatrixMarketFile("Data/mask.mm", *Wp);
if(err!=0){ 
  MPI_Finalize();
  exit(1);
}

    //(*Wp).ReplaceRowMap(Wp_Map);
/*    cout << "size of pmask is " << (*Wp).NumGlobalRows() << " by " << (*Wp).NumGlobalCols() << ", numzeros: " << (*Wp).NumGlobalNonzeros() << endl;
    cout << "size of pmask is " << (*Wp).NumMyRows() << " by " << (*Wp).NumMyCols() << ", numzeros: " << (*Wp).NumMyNonzeros() << endl;
*/
    err = ConstructWFOperator(n, A, frm_comm, Wp);

    //cout << "size of masked operator is " << (*A).NumGlobalRows() << " by " << (*A).NumGlobalCols() << ", nonzeros is " << (*A).NumGlobalNonzeros()  << endl;




    // define the map for the reading in data and load the data
    Epetra_Map B_range(n*n, n*n, 0, frm_comm);
    const char* matFile = "grads16.128.mat";


    //each proc gets the name of its mat variables, based on frame
    std::stringstream sstmx;
    sstmx << "x" << MyFrame;
    string xvar = sstmx.str();

    std::stringstream sstmy;
    sstmy << "y" << MyFrame;
    string yvar = sstmy.str();

    char *varNameX;
    char *varNameY;
    if(MyFrame < 10){
    varNameX = new char[3];
    varNameY = new char[3];
    strncpy(varNameX, xvar.c_str(), 3 );
    strncpy(varNameY, yvar.c_str(), 3 );
    } else{ //handles up to 99 frames
    varNameX = new char[4];
    varNameY = new char[4];
    strncpy(varNameX, xvar.c_str(), 4 );
    strncpy(varNameY, yvar.c_str(), 4 );
    }
//cout << "Process " << MyPID << " reporting, loading " << varNameX << " and " << varNameY << endl;
    //err = readMat(matFile, varName, xGradMeas, yGradMeas, range);
    comm.Barrier();
    Epetra_Time timeSolver(frm_comm);
    err = readMat(matFile, varNameX, xGradMeas, B_range);
if(err!=0){
MPI_Finalize(); exit(1);}
    err = readMat(matFile, varNameY, yGradMeas, B_range);
if(err!=0){
MPI_Finalize(); exit(1);}


    // BELOS
    // Some typedefs for oft-used data types
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;
    typedef Belos::MultiVecTraits<double, Epetra_MultiVector> MVT;
    bool ierr;

    // Define new operator.
    (*A).OptimizeStorage();
    Teuchos::RCP<OP> MyOp(*&A);
        
    if (verbose) {
      cout << "Belos Example: LSQR " << endl;
      cout << "Problem info:" << endl;
      cout << endl;
     }
    
    // Setup some more problem/solver parameters:
    int blocksize = 1;
  
    const Epetra_Map * XMap = &(A->DomainMap());
    Teuchos::RCP<Epetra_MultiVector> Bx = Teuchos::rcp( new Epetra_MultiVector(*xGradMeas));
    Teuchos::RCP<Epetra_MultiVector> By = Teuchos::rcp( new Epetra_MultiVector(*yGradMeas));
//cout << (*A).NormFrobenius() << " is norm of A" << endl;

    // Stack Bx and By
    Epetra_Map BMap(2*n*n, 0, frm_comm);
    Teuchos::RCP<Epetra_MultiVector> B = Teuchos::rcp( new Epetra_MultiVector(BMap,blocksize) );
//    cout << "local len, global len:" << endl;
//    cout << B->MyLength() << ", " << B->GlobalLength() << endl;
    int n2 = n*n;
    double * bRes;

    //Stack x- and y-gradient vectors into the B vector
    int nZB;
    double *BxVals = (double*) malloc(n2*sizeof(double));
    double *ByVals = (double*) malloc(n2*sizeof(double));

    BxVals = (*Bx)[0];
    ByVals = (*By)[0];

    int *BxIdx = new int[n2];
    int *ByIdx = new int[n2];
    for(i=0; i < n2; i++){
      B->ReplaceGlobalValue(i     , 0, BxVals[i]);
      B->ReplaceGlobalValue(i + n2, 0, ByVals[i]);
      //assert((*B)[0][i] == (*xGradMeas)[0][i]);
      //assert((*B)[0][i + n2] == (*yGradMeas)[0][i]);
      
    }
//MultiVectorToMatrixMarketFile("Data/B.mm", *B);

    Teuchos::RCP<Epetra_MultiVector> X = Teuchos::rcp( new Epetra_MultiVector(*XMap,blocksize) );
    // Initialize the solution with zero and right-hand side with random entries
    X->PutScalar( 0.0 );
            //err = MultiVectorToMatrixMarketFile(  "Data/B.mm", *B );
            //err = RowMatrixToMatrixMarketFile(  "Data/Aop.mm", *A );
    
            //err = MultiVectorToMatrixMarketFile(  "Data/X.mm", *X );
//            err = MultiVectorToMatrixMarketFile(  "Data/Bx.mm", *Bx );
//            err = MultiVectorToMatrixMarketFile(  "Data/By.mm", *By );
    // Setup the linear problem, with the matrix A and the vectors X and B
    Teuchos::RCP< Belos::LinearProblem<double,MV,OP> > myProblem = Teuchos::rcp( new Belos::LinearProblem<double,MV,OP>(MyOp, X, B) );
    myProblem->setLeftPrec(Teuchos::null);
 
    // Signal that we are done setting up the linear problem
    ierr = myProblem->setProblem();

    // Check the return from setProblem(). 
    assert(ierr == true);

    // Specify the verbosity level. Options include:
    // Belos::Errors 
    //   This option is always set
    // Belos::Warnings 
    //   Warnings (less severe than errors)
    //
    // Belos::IterationDetails 
    //   Details at each iteration, such as the current eigenvalues
    // Belos::OrthoDetails 
    //   Details about orthogonality
    // Belos::TimingDetails
    //   A summary of the timing info for the solve() routine
    // Belos::FinalSummary 
    //   A final summary 
    // Belos::Debug 
    //   Debugging information
    int verbosity = Belos::Warnings + Belos::Errors + Belos::FinalSummary + Belos::TimingDetails;
    int lessverbose = Belos::Errors + Belos::TimingDetails; 
    int none = 0;
    // Create the parameter list for the eigensolver
    Teuchos::RCP<Teuchos::ParameterList> myPL = Teuchos::rcp( new Teuchos::ParameterList() );
    //myPL->set( "Verbosity", verbosity);
    myPL->set( "Verbosity", lessverbose);
    myPL->set( "Block Size", blocksize );
    myPL->set( "Maximum Iterations",500 );
    myPL->set( "Lambda", 0.001); // regularization parameter
    // myPL->set( "Rel RHS Err", 0.001 ); // stopping criterior parameter, default = 0.
    myPL->set( "Rel Mat Err", 0.001 ); // stopping criterior parameter, equivalent to 'tol' in MATLAB lsqr.
  

    // Create the LSQR solver
    // This takes as inputs the linear problem and the solver parameters
    Belos::LSQRSolMgr<double, MV, OP> mySolver(myProblem, myPL);
    Teuchos::RCP<const Teuchos::ParameterList> valParams = mySolver.getValidParameters();
//    cout << "Valid Params: " << endl;
//    valParams->print(cout, 0, 1, 1);

    // Solve the linear problem, and save the return code
    Belos::ReturnType solverRet = mySolver.solve();
    comm.Barrier();
    double timing = timeSolver.ElapsedTime();

    // Check return code of the solver: Unconverged, Failed, or OK
    switch (solverRet) {
    // UNCONVERGED
        case Belos::Unconverged:
        if (lead_proc)
            cout << "Belos::BlockLSQRSolMgr::solve() did not converge!" << endl;
    // output the computed solution.
        err = MultiVectorToMatrixMarketFile(  "Data/unconv_sol.mm", *X );
        
        #ifdef HAVE_MPI
         MPI_Finalize();
        #endif
        return 0;
        break;

     // CONVERGED
     case Belos::Converged:
        if (lead_proc){
            cout << "Belos::BlockLSQRSolMgr::solve() converged!" << endl;
            std::stringstream outstm;
            outstm << "Data/conv_sol_frm" << MyFrame << ".mm" ;
            string outfile = outstm.str();
            err = MultiVectorToMatrixMarketFile(  outfile.c_str(), *X );
        }
        break;
    }
  
    double iterResNorm = mySolver.getResNorm();
    std::vector<double> normB(blocksize);
    MVT::MvNorm( *Bx, normB );

    // Output results to screen
    if(verbose) {
    //cout << scientific << setprecision(6) << showpoint;
        cout << "******************************************************\n"
             << "           Results (outside of linear solver)           \n"
             << "------------------------------------------------------\n"
             << "  Linear System\t\tRelative Residual\n"
             << "------------------------------------------------------\n";
        for( int i=0 ; i<blocksize ; ++i ) {
          cout << "  " << i+1 << "\t\t\t" << iterResNorm/normB[i] << endl;
        }

        cout << "******************************************************\n" << endl;
    }
    comm.Barrier();
    if(lead_proc){
        cout << "Time to solve frame " << MyFrame << ": " << timing << endl; 
        //Teuchos::Array<Teuchos::RCP<Teuchos::Time> > belosTime = myProblem->getTimers();
    }

    /*    if(verbose)
        cout << "the computed residul norm is " << mySolver.getResNorm() << endl;

    // Test residuals
    Epetra_MultiVector Rresidual(xGradMeas->Map(), blocksize );
    MyOperator->Apply( *X, Rresidual ); // R = A*X
    MVT::MvAddMv( -1.0, *Bx, 1.0, Rresidual, Rresidual );// R -= B 

    // Compute the 2-norm of each vector in the MultiVector
    // and store them to a std::vector<double>
    std::vector<double> normR(blocksize), normB(blocksize);
    MVT::MvNorm( Rresidual, normR );
    if(verbose)
        cout << "  " << normR[0]/normB[0] << endl; */
    
//    delete A; //causes free() to segfault upon exit
//    delete Wp;
    delete xGradMeas;
    delete yGradMeas;

    #ifdef HAVE_MPI
        MPI_Finalize();
    #endif

    

    return(EXIT_SUCCESS);
}
