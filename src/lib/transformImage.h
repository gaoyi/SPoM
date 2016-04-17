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


#ifndef transformImage_h_
#define transformImage_h_


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// itk
#include "itkVersorRigid3DTransform.h"
#include "itkSimilarity3DTransform.h"


namespace ShapeAnalysis
{
  /**
   * SIMILARITY Transform the image. 
   * 
   * interpolationType = 0 for NN interp, 1 for linear interp, 2 for
   * bspline interp. Default is 1
   */
  template<typename InputImageType, typename ReferenceImageType, typename OutputImageType>
  typename OutputImageType::Pointer
  transformImage(typename itk::Similarity3DTransform<double>::Pointer transform, \
                 typename InputImageType::Pointer inputImage, typename ReferenceImageType::Pointer referenceImage, \
                 typename OutputImageType::PixelType fillInValue, char interpolationType = 1);


  /**
   * RIGID Transform the image. 
   * 
   * interpolationType = 0 for NN interp, 1 for linear interp, 2 for
   * bspline interp. Default is 1
   */
  template<typename InputImageType, typename ReferenceImageType, typename OutputImageType>
  typename OutputImageType::Pointer
  transformImage(typename itk::VersorRigid3DTransform<double>::Pointer transform, \
                 typename InputImageType::Pointer inputImage, typename ReferenceImageType::Pointer referenceImage, \
                 typename OutputImageType::PixelType fillInValue, char interpolationType = 1);


  /**
   * Warp image using vector image
   * 
   * interpolationType = 0 for NN interp, 1 for linear interp, 2 for
   * bspline interp.
   */
  template<typename ReferenceImageType, typename MovingImageType, typename DisplacementFieldType>
  typename ReferenceImageType::Pointer
  warpImage(typename ReferenceImageType::Pointer img,                     \
            typename MovingImageType::Pointer movingImg,                  \
            typename DisplacementFieldType::Pointer displacementField, \
            typename ReferenceImageType::PixelType fillInVal, \
            char interpolationType);


}// namespace ShapeAnalysis

#include "transformImage.hxx"


#endif //transformImage_h_
