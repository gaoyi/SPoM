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


#ifndef PCGSPDLinearEquationSolver_h_
#define PCGSPDLinearEquationSolver_h_

// itk
#include "vnl/vnl_vector.h"
#include "vnl/vnl_sparse_matrix.h"

namespace ShapeAnalysis
{
  // Solvie Ax = b using Preconditioned Conjugate Gradient method.
  // A MSUT be Symetric-positive-definite and sparse, though not
  // checked.
  template<typename DataType>
  void
  PCGSPDLinearEquationSolver(const vnl_sparse_matrix<DataType>& A, const vnl_vector<DataType>& b, \
                             vnl_vector<DataType>& x, \
                             DataType tol = 1e-10);
} //


#include "PCGSPDLinearEquationSolver.hxx"

#endif
