/*******************************************************************************/
/*                                                                             */
/*  Copyright (c) 2010, Yi Gao                                                 */
/*  gaoyi@gatech.edu                                                           */
/*                                                                             */
/*                                                                             */
/*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR */
/*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   */
/*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    */
/*  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER */
/*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    */
/*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        */
/*  DEALINGS IN THE SOFTWARE.                                                  */
/*                                                                             */
/*  See the README.md and COPYING files for usage and copyright information.   */
/*                                                                             */
/*******************************************************************************/


#ifndef PCGSPDLinearEquationSolver_hxx_
#define PCGSPDLinearEquationSolver_hxx_

#include <iostream>

// itk
#include "vnl/vnl_vector.h"
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_math.h"

// local
#include "PCGSPDLinearEquationSolver.h"

namespace ShapeAnalysis
{
  template<typename DataType>
  void
  PCGSPDLinearEquationSolver(const vnl_sparse_matrix<DataType>& A, const vnl_vector<DataType>& b, \
                             vnl_vector<DataType>& x, \
                             DataType tol)
  {
    // Solvie Ax = b using Preconditioned Conjugate Gradient method.
    // A MSUT be Symetric-positive-definite and sparse, though not
    // checked.
    typedef vnl_vector<DataType> VNLVectorType;

    unsigned long numUnknowns = A.cols();

    if (x.size() != numUnknowns)
      {
        x.set_size(numUnknowns);
        x.fill(0.0);
      }

    VNLVectorType rx = b;
    VNLVectorType zx(numUnknowns);

    // Jacobi preconditioner
    VNLVectorType Dinv(numUnknowns);
    for (unsigned long ip = 0; ip < numUnknowns; ++ip )
      {
        Dinv[ip] = 1.0 / ( A(ip, ip) + vnl_math::eps );

        zx[ip] = rx[ip] * Dinv[ip];
      }

    VNLVectorType dx = zx;

    long numIter = b.size();
    if ( b.size() != numUnknowns )
      {
        // check for safe
        std::cerr << "b.size() != numUnknowns\n";
      }
    numIter += numIter / 10; 
    // let the iteration times a little more than the dimension

    //double tol = 1e-12;

//    std::cout<<"in iteration\n";
    for ( long i = 0; i <= numIter; ++i )
      {
        //std::cout<<i<<", "<<std::flush;

        VNLVectorType Dxd;
        A.pre_mult(dx, Dxd);

        double dDxd = inner_product(dx, Dxd);

        double zxTrx = inner_product(zx, rx);

        double alphax = zxTrx / ( dDxd + vnl_math::eps );

        x += (dx*alphax);

        rx -= (Dxd*alphax);

        double rxTrx = inner_product(rx, rx);
        if ( rxTrx < tol)
          {
//            std::cout<<"out from here when i = "<<i<<std::endl;
            break;
          }

        for ( unsigned long id = 0; id < numUnknowns; ++id )
          {
            zx[id] = rx[id] * Dinv[id];
          }

        double betaX = inner_product(zx, rx) / ( zxTrx + vnl_math::eps );

        dx = zx + (dx*betaX);
      }
//    std::cout<<std::endl;

    return;
  }
} //namespace ShapeAnalysis


#endif
