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


#ifndef SumOfTwoPoissonShapeAnalysis_h_
#define SumOfTwoPoissonShapeAnalysis_h_

#include <vector>

// itk
#include "itkImage.h"

// vtk
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace ShapeAnalysis
{
  template< typename TShapeImageType >
  class SumOfTwoPoissonShapeAnalysis
  {
  public:
    static const unsigned int Dim = TShapeImageType::ImageDimension;

    typedef TShapeImageType ShapeImageType;
    typedef SumOfTwoPoissonShapeAnalysis< ShapeImageType > Self;

    typedef typename ShapeImageType::Pointer ShapeImageTypePointer;

    typedef double FloatType;
    typedef itk::Image<FloatType, Dim> FloatImageType;
    typedef typename FloatImageType::Pointer FloatImagePointer;

    SumOfTwoPoissonShapeAnalysis();
    ~SumOfTwoPoissonShapeAnalysis() {}

    void SetShapeGroup1(const std::vector< ShapeImageTypePointer >& shapeGroup);
    void SetShapeGroup2(const std::vector< ShapeImageTypePointer >& shapeGroup);

    void setMeanShapeImage(FloatImagePointer meanshape);

    FloatImagePointer getMeanShapeImage();
    vtkSmartPointer<vtkPolyData> getMeanShapeSurface();

    void update();

    void setPerformRegistration(bool reg);

  private:
    // data
    std::vector< ShapeImageTypePointer > m_shapeGroup1;
    std::vector< ShapeImageTypePointer > m_shapeGroup2;

    std::vector< FloatImagePointer > m_similarityRegisteredShapeGroup1;
    std::vector< FloatImagePointer > m_similarityRegisteredShapeGroup2;

    std::vector< FloatImagePointer > m_poissonMapGroup1;
    std::vector< FloatImagePointer > m_poissonMapGroup2;

    FloatImagePointer m_registrationFixedImage;

    FloatImagePointer m_meanShape;
    vtkSmartPointer<vtkPolyData> m_meanShapeSurface;

    std::vector< std::vector<FloatType> > m_poissonDistanceOnMeanSurface1;
    std::vector< std::vector<FloatType> > m_poissonDistanceOnMeanSurface2;

    std::vector<FloatType> m_pValuesOnMeanShapeSurface;
    std::vector<FloatType> m_pValuesOnMeanShapeSurfaceFDR;

    bool m_performRegistration;

    bool m_useOutsideMeanShape;

    bool m_allDone;

    // fn
    void _preprocess();
    void _cropInputImages();
    void _isotropicizeFirstImage();
    void _computeMeanShape();
    void _computeMeanShapeSurface();
    void _computePoissonMaps();
    void _poissonDistanceMapOnMeanShapeSurface();
    void _poissonDistance(FloatImagePointer poissonField, double startingPosition[3]); // poissonField should be const pointer indeed
    void _computePValueOnMeanShape();
  };

}// ShapeAnalysis


#include "SumOfTwoPoissonShapeAnalysis.hxx"

#endif
