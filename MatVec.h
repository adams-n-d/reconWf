#ifndef _AL_SPMV_H
#define _AL_SPMV_H

#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>

#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include "Epetra_Operator.h"

#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_BlockMapIn.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_VectorIn.h>

#include "Epetra_Util.h"
#include "Epetra_BlockMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "EpetraExt_mmio.h"

#include "BelosConfigDefs.hpp"
#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosTypes.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosLSQRSolMgr.hpp"

class SPMV: public virtual Epetra_Operator {
private:
  typedef Epetra_MultiVector                MV;
  typedef Epetra_CrsMatrix                MTRX;
  typedef Epetra_Map                       MAP;

    MTRX A;
    MTRX Wp;
    MTRX R;
    MAP  *mapA;
    MAP  *mapWp;
    MAP  *mapR;
 
  bool UseTranspose_;

public:

  SPMV( const MTRX &A_, const MTRX &Wp_, const MTRX &R_):
    A(A_), Wp(Wp_), R(R_){
        mapA = new MAP(A.DomainMap());
        mapWp = new MAP(Wp.RangeMap());
        mapR = new MAP(R.RangeMap());
    }

  ~SPMV(){};

  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);}

  int Apply(const Epetra_MultiVector &, Epetra_MultiVector &) const;

  int ApplyInverse( const Epetra_MultiVector &, Epetra_MultiVector &) const{
    return (-1); // Not implemented.
  }

  double NormInf() const {
    return(-1);
  }

  const char* Label() const {
    return("Q verstion of matrix-vector multiplication");
  }

  bool UseTranspose() const {return(UseTranspose_);}

  bool HasNormInf () const {
    return(false);
  }

  const Epetra_Comm & Comm() const {
    return(A.Comm());
  }

  const Epetra_Map & OperatorDomainMap() const {
     // This the sentence which bother me most when print out the output info.
    //cout << "Domain map" << endl << *mapA << endl;
    return(*mapA);
  }

  const Epetra_Map & OperatorRangeMap() const {
     // This the sentence which bother me most when print out the output info.
    //cout << "Range map" << endl << *mapA << endl;
    return(*mapR);
  }

  // Other function.
  int Multiply( bool , const Epetra_MultiVector &, Epetra_MultiVector &) const;

};

/*  int SPMV::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Multiply(UseTranspose_, X, Y));}
	
  int SPMV::Multiply( bool TransA, const Epetra_MultiVector &x, Epetra_MultiVector &y) const {

    int err;
    //int k, l;
    int numxVec = x.NumVectors();

    if(!TransA) {

        // H*x.
        MV Hx(MV(H.RangeMap(), numxVec));
        err = H.Multiply(false, x, Hx);
        if (err != 0)
            return err;

        // A*H*x = A*temp.
        MV AHx(MV(A.RangeMap(), numxVec));
        err = A.Multiply(false, Hx, AHx);
        if (err != 0)
            return err;

        y = AHx;
    }
    else {

        // ATx.
        MV ATx(MV(A.DomainMap(), numxVec));
        err = A.Multiply(true, x, ATx);
        if (err != 0)
            return err;

        //HT*AT*x = HT*temp
        MV HTATx(MV(H.DomainMap(), numxVec));
        err = H.Multiply(true, ATx, HTATx);
        if (err != 0)
            return err;

        y = HTATx;
    }

    return 0;
}*/
#endif
