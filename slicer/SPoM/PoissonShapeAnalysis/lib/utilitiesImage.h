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


#ifndef utilitiesImage_h_
#define utilitiesImage_h_

// itk
#include "itkAffineTransform.h"
#include "itkImage.h"

// vnl
#include "vnl/vnl_matrix.h"


namespace ShapeAnalysis
{
  /**
   * Cast image pixel type
   */
  template< typename inputPixel_t, typename outputPixel_t > 
  typename itk::Image<outputPixel_t, 3 >::Pointer
  castItkImage( typename itk::Image<inputPixel_t, 3>::Pointer inputImage );

  template<typename image_t>
  double getVol(typename image_t::Pointer img, typename image_t::PixelType thld = 0);

  /**
   * binarilize image
   */
  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                 \
                  typename input_image_t::PixelType thld,                         \
                  typename output_image_t::PixelType insideValue = 1);


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                 \
                  typename input_image_t::PixelType lowerT,                         \
                  typename input_image_t::PixelType upperT, \
                  typename output_image_t::PixelType insideValue,               \
                  typename output_image_t::PixelType outsideValue);


  template<typename ImageType>
  void binarilizeImageSeries(std::vector<typename ImageType::Pointer>& imgList);



  /**
   * Compute the non-zero region of the image
   */
  template<typename image_t>
  typename image_t::RegionType
  computeNonZeroRegion(typename image_t::Pointer img);

  /**
   * Enlarge the non-zero region so that the region is not too tightly around the non-zero reigon
   */
  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegion(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion);


  /**
   * Crop the mask by its non-zero region
   */
  template<typename MaskImageType >
  typename MaskImageType::Pointer
  cropNonZeroRegionFromImage(typename MaskImageType::Pointer mask);


  /*********************************************************************************
   * About atals segmentation
   *********************************************************************************/

  /**
   * Extract the ROI from the image using the region
   */
  template<typename image_t>
  typename image_t::Pointer
  extractROI(typename image_t::Pointer img, typename image_t::RegionType region);


  /** 
   * Generate an all-zero image the same size/origin/spacing/etc. as
   * referenceImg, inside of whick, the roiRegion is the roiImg
   */
  template<typename image_t>
  typename image_t::Pointer
  antiExtractROI(typename image_t::ConstPointer roiImg,         \
                 const typename image_t::RegionType roiRegion,  \
                 typename image_t::ConstPointer referenceImg);


  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegionByOnePixel(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion);

  /**
   * convert between probability image and logit image 
   */
  template<typename TInputImage, typename TOutputImage>
  typename TOutputImage::Pointer
  probabilityImageToLogitImage(typename TInputImage::Pointer inputImage);

  template<typename TInputImage, typename TOutputImage>
  typename TOutputImage::Pointer
  logitImageToProbabilityImage(typename TInputImage::Pointer inputImage);

}// namespace ShapeAnalysis


#include "utilitiesImage.hxx"

#endif
