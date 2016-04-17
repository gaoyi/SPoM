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


#ifndef ConformalMetricNeumannDirichletPoissonSolver3D_h_
#define ConformalMetricNeumannDirichletPoissonSolver3D_h_

#include <vector>

// itk
#include "itkImage.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_sparse_matrix.h"

namespace ShapeAnalysis
{

  /**
   * Solve the Poisson eq with both Neumann and Dirichelet BC. The
   * input are two images indicating the Dirichlet BC location and the
   * values. Specifically, the Non-zero voxel of the DBC-mask image
   * means that the BC at that location is the value of the BC-value
   * image. For the boundary, we use the Neumann BC. This is the
   * difference between this and the DirichletPoissonSolver3D class.
   */
  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  class ConformalMetricNeumannDirichletPoissonSolver3D
  {
  public:
    static const unsigned int Dim = TNeumannDirichletBoundaryConditionMaskType::ImageDimension;

    typedef TInputImageType InputImageType;

    typedef TNeumannDirichletBoundaryConditionMaskType DirichletBoundaryConditionMaskType;
    typedef DirichletBoundaryConditionMaskType MaskImageType;

    typedef TNeumannDirichletBoundaryConditionValueImageType DirichletBoundaryConditionValueImageType;
    typedef DirichletBoundaryConditionValueImageType BCValueImageType;

    typedef ConformalMetricNeumannDirichletPoissonSolver3D< InputImageType, MaskImageType, BCValueImageType > Self;

    typedef typename MaskImageType::IndexType IndexType;

    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename BCValueImageType::Pointer BCValueImagePointer;

    typedef double FloatType;
    typedef itk::Image<FloatType, Dim> FloatImageType;
    typedef typename FloatImageType::Pointer FloatImagePointer;

    ConformalMetricNeumannDirichletPoissonSolver3D();
    ~ConformalMetricNeumannDirichletPoissonSolver3D() {}

    void setInputImage(const InputImageType* inputImage);

    void setBoundaryConditionMask(MaskImagePointer bcMask);
    void setBoundaryConditionValue(BCValueImagePointer bcValueImage);

    void setBeta(double b);

    void update();

    FloatImagePointer getPoissonImage();

  private:
    // data
    const InputImageType* m_inputImage; // the conformal metric is derived from the gradient of this image.

    MaskImagePointer m_DBCMask;
    BCValueImagePointer m_DBCValueImage; //TODO: for random walker, this will be drived from initial label map, so this can be computed and not need to be given. While, another idea is to keep this Possion solver as general "conformal Poisson solver" and before this there is an RW class just do the pre-processing and then call this. Which one is better?

    double m_beta;

    FloatImagePointer m_poissonImage;

    const FloatImageType* m_initialPossionImage;

    long m_numberOfDBCVoxel;
    long m_numberOfUnknownVoxel;

    typedef vnl_vector<FloatType> VNLVectorType;

    VNLVectorType m_XKnown;
    VNLVectorType m_XUnknown;
    vnl_sparse_matrix<FloatType> m_B;
    vnl_sparse_matrix<FloatType> m_S;

    typedef itk::Image<long, Dim> LongImageType;
    typename LongImageType::Pointer m_indexMap;
    // for each voxel, what's its index in the long vector? the BC
    // voxels are at the front of the long vector, followed by the
    // unknowns.

    bool m_allDone;

    // fn
    void _preprocess();
    void _buildIndexMapping();
    void _setupEquation();
    void _solve();
  };

}// ShapeAnalysis


#include "ConformalMetricNeumannDirichletPoissonSolver3D.hxx"

#endif
