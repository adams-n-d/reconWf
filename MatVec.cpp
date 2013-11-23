#include "MatVec.h"

int SPMV::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Multiply(UseTranspose_, X, Y));
}

int SPMV::Multiply( bool TransA, const Epetra_MultiVector &x, Epetra_MultiVector &y) const {

    int err;
    //int k, l;
    int numVec = x.NumVectors();

    if(!TransA) {

        // A*x.
        //numVec = A.NumGlobalRows();
        MV Ax(MV(A.RangeMap(), numVec));
        err = A.Multiply(false, x, Ax);
        if (err != 0)
            return err;

        // Wp*A*x = Wp*temp.
        //numVec = Wp.NumGlobalRows();
        MV WAx(MV(Wp.RangeMap(), numVec));
        err = Wp.Multiply(false, Ax, WAx);
        if (err != 0)
            return err;

        // R*Wp*A*x.
        //numVec = R.NumGlobalRows();
        MV RWAx(MV(R.RangeMap(),numVec));
        err = R.Multiply(false,WAx,RWAx);
        if(err!=0)
            return err;

        y = RWAx;
    }
    else {

        // RTx.
        //numVec = R.NumGlobalCols(); 
        MV RTx(MV(R.DomainMap(), numVec));
        err = R.Multiply(true, x, RTx);
        if (err != 0)
            return err;

        //WT*RT*x = WT*temp
        //numVec = Wp.NumGlobalCols();
        MV WTRTx(MV(Wp.DomainMap(), numVec));
        err = Wp.Multiply(true, RTx, WTRTx);
        if (err != 0)
            return err;

        //AT*WT*RT*x.
        //numVec = A.NumGlobalCols();
        MV ATWTRTx(MV(A.DomainMap(),numVec));
        err = A.Multiply(true,WTRTx,ATWTRTx);
        if(err!=0)
            return err;

        y = ATWTRTx;
    }

    return 0;
}
