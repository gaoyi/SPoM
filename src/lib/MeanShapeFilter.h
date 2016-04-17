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


#ifndef MeanShapeFilter_h_
#define MeanShapeFilter_h_

#include <vector>

// itk
#include "itkImage.h"

namespace ShapeAnalysis
{

  /**
   * Given a set of input shapes, compute their mean shape. The input
   * shape can be given as floating point element type. A threshold is
   * needed to get the binary from the input shapes. Two ways is used
   * to get the mean shape: average of the binary mean shapes, and the
   * log-odd method. The output volumes are 0-1 binary, regardless of
   * the pixel type.
   *
   * The input shapes MUST have the same regions.
   */
  template< typename TInputShapeType, typename TOutputShapeType >
  class MeanShapeFilter
  {
  public:
    static const unsigned int Dim = TInputShapeType::ImageDimension;

    typedef TInputShapeType InputShapeType;
    typedef TOutputShapeType OutputShapeType;

    typedef MeanShapeFilter< InputShapeType, OutputShapeType > Self;

    typedef typename InputShapeType::IndexType IndexType;
    typedef typename InputShapeType::RegionType RegionType;

    typedef typename InputShapeType::Pointer InputShapePointer;
    typedef typename OutputShapeType::Pointer OutputShapePointer;

    typedef typename InputShapeType::PixelType InputPixelType;
    typedef typename OutputShapeType::PixelType OutputPixelType;

    typedef double FloatType;
    typedef itk::Image<FloatType, Dim> FloatImageType;
    typedef typename FloatImageType::Pointer FloatImagePointer;

    MeanShapeFilter();
    ~MeanShapeFilter() {}

    void setInputShapes(const std::vector<InputShapePointer>& inputShapes);
    void setThresholds(FloatType thldLower, FloatType thldUpper);
    void setMeanShapeType(char meanShapeType);

    void update();

    FloatImagePointer getFloatMeanShape();
    OutputShapePointer getMeanShape();

  private:
    // data
    std::vector<InputShapePointer> m_inputShapes;

    FloatType m_inputThresholdUpper;
    FloatType m_inputThresholdLower;

    FloatType m_epsilonForHeviside;

    FloatImagePointer m_floatMeanShape;
    OutputShapePointer m_meanShape;

    /**
     * 0: log-odds mean shape (default)
     * 1: average of binaries  
     */
    char m_meanShapeType;

    bool m_allDone;

    // fn
    void _computeAverageBinaryMeanShape();
    void _computeLogoddsMeanShape();
  };

}// ShapeAnalysis


#include "MeanShapeFilter.hxx"

#endif
