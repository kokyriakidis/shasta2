#include "SHASTA2_ASSERT.hpp"
#include "shastaLapack.hpp"
using namespace shasta2;

#include "algorithm.hpp"



// dgesvd wrapper for boost ublas matrix.
void shasta2::dgesvd(
    boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& A,
    vector<double>& S,
    boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& U,
    boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& VT
    )
{

    char JOBU = 'A';
    char JOBVT = 'A';

    int M = int(A.size1());
    int N = int(A.size2());
    SHASTA2_ASSERT(M > 0);
    SHASTA2_ASSERT(N > 0);

    int LDA = M;

    S.resize(min(M, N));

    U.resize(M, M);
    int LDU = M;

    VT.resize(N, N);
    int LDVT = N;

    int LWORK = 10 * max(M, N);
    vector<double> WORK(LWORK);

    int INFO = 0;

    dgesvd_(&JOBU, &JOBVT, M, N, &A(0, 0), LDA, &S[0], &U(0, 0), LDU, &VT(0, 0), LDVT, &WORK[0], LWORK, INFO);
    SHASTA2_ASSERT(INFO == 0);
}
