#ifndef SHASTA_SHASTA_LAPACK_HPP
#define SHASTA_SHASTA_LAPACK_HPP


/*******************************************************************************

Shasta uses some functionality provided by the Lapack and Blas libraries.

Lapack and Blas come with this license:

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
 .
 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 .
 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer listed
   in this license in the documentation and/or other materials
   provided with the distribution.
 .
 - Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
 .
 The copyright holders provide no reassurances that the source code
 provided does not infringe any patent, copyright, or any other
 intellectual property rights of third parties.  The copyright holders
 disclaim any liability to any recipient for claims brought against
 recipient by any third party for infringement of that parties
 intellectual property rights.
 .
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*******************************************************************************/

#include <boost/numeric/ublas/matrix.hpp>

#include <vector.hpp>

extern "C" void dgemv_(
    const char* TRANS,
    const int& M,
    const int& N,
    const double& ALPHA,
    const double* A,
    const int& LDA,
    double* X,
    const int& INCX,
    const double& BETA,
    const double* Y,
    const int& INCY
    );

extern "C" void dgesvd_(
    const char* JOBU,
    const char* JOBVT,
    const int& M,
    const int& N,
    double* A,
    const int& LDA,
    double* S,
    double* U,
    const int& LDU,
    double* VT,
    const int& LDVT,
    double* WORK,
    const int& LWORK,
    int& INFO
    );

extern "C" void dsyev_(
    const char* JOBZ,
    const char* UPLO,
    const int& N,
    double* A,
    const int& LDA,
    double* W,
    double* WORK,
    const int& LWORK,
    int& INFO
    );


namespace shasta {

    // dgesvd wrapper for boost ublas matrix.
    void dgesvd(
        boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& A,
        vector<double>& S,
        boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& U,
        boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>& VT
        );
}

#endif
