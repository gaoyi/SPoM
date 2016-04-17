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


#ifndef ShapeSumOfTwoPoissonFilter_h_
#define ShapeSumOfTwoPoissonFilter_h_

#include <vector>

// itk
#include "itkImage.h"

// vtk
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace ShapeAnalysis
{
  /**
   * The inpue MUST be 0-1 integer typed image.
   */
  template< typename TShapeImageType >
  class ShapeSumOfTwoPoissonFilter
  {
  public:
    static const unsigned int Dim = TShapeImageType::ImageDimension;

    typedef TShapeImageType ShapeImageType;
    typedef ShapeSumOfTwoPoissonFilter< ShapeImageType > Self;

    typedef typename ShapeImageType::Pointer ShapeImagePointer;

    typedef double FloatType;
    typedef itk::Image<FloatType, Dim> FloatImageType;
    typedef typename FloatImageType::Pointer FloatImagePointer;

    ShapeSumOfTwoPoissonFilter();
    ~ShapeSumOfTwoPoissonFilter() {}

    void setInputShape(ShapeImagePointer inputShape);
    //void setNumberOfIterations(long n);

    void update();

    FloatImagePointer getPoissonImage();

  private:
    // data
    ShapeImagePointer m_inputShape;

    ShapeImagePointer m_boundaryMask;
    ShapeImagePointer m_kernelMask;

    typedef itk::Image<char, Dim> CharImageType;
    typename CharImageType::Pointer m_DirichletBCMaskForPoissonSolver;

    FloatImagePointer m_DirichletBCForPoissonSolver;

    FloatImagePointer m_poissonImage;

    FloatImagePointer m_dfdt;

    FloatType m_dt;

    //long m_numIter;

    bool m_allDone;

    // fn
    //void _preprocess();
    void _computeInsidePoisson();
    void _computeOutsidePoisson();
    void _computeDirichletBC();
    void _solvePoisson();
    void _computeSumOfTwoPoisson();

    long _numberOfNonZeroVoxels(ShapeImagePointer img);

  };

}// ShapeAnalysis


#include "ShapeSumOfTwoPoissonFilter.hxx"

#endif
