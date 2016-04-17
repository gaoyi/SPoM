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


#ifndef DirichletPoissonSolver3D_h_
#define DirichletPoissonSolver3D_h_

#include <vector>

// itk
#include "itkImage.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_sparse_matrix.h"

namespace ShapeAnalysis
{

  /**
   * Solve the Poisson eq with Dirichelet BC. The input is a single
   * image indicating the BC. Specifically, the Non-zero voxel of the
   * BC-mask image means that the BC at that location is the value of the
   * BC-value image. Moreover, the boundary is forced to have 0 value, unless
   * they have other values in the BC image. In that case, the BC
   * image boundary values will be used.
   */
  template< typename TDirichletBoundaryConditionMaskType, typename TDirichletBoundaryConditionValueImageType >
  class DirichletPoissonSolver3D
  {
  public:
    static const unsigned int Dim = TDirichletBoundaryConditionMaskType::ImageDimension;

    typedef TDirichletBoundaryConditionMaskType DirichletBoundaryConditionMaskType;
    typedef DirichletBoundaryConditionMaskType MaskImageType;

    typedef TDirichletBoundaryConditionValueImageType DirichletBoundaryConditionValueImageType;
    typedef DirichletBoundaryConditionValueImageType BCValueImageType;

    typedef DirichletPoissonSolver3D< MaskImageType, BCValueImageType > Self;

    typedef typename MaskImageType::IndexType IndexType;

    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename BCValueImageType::Pointer BCValueImagePointer;

    typedef double FloatType;
    typedef itk::Image<FloatType, Dim> FloatImageType;
    typedef typename FloatImageType::Pointer FloatImagePointer;

    DirichletPoissonSolver3D();
    ~DirichletPoissonSolver3D() {}

    void setBoundaryConditionMask(MaskImagePointer bcMask);
    void setBoundaryConditionValue(BCValueImagePointer bcValueImage);

    void update();

    FloatImagePointer getPoissonImage();

  private:
    // data
    MaskImagePointer m_BCMask;
    BCValueImagePointer m_BCValueImage;

    FloatImagePointer m_poissonImage;

    long m_numberOfBCVoxel;
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


#include "DirichletPoissonSolver3D.hxx"

#endif
