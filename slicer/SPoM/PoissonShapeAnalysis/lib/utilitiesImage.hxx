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


#ifndef utilitiesImage_hxx_
#define utilitiesImage_hxx_

#include <csignal>

// itk
#include "itkAffineTransform.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkDanielssonDistanceMapImageFilter.h"

#include "itkImageDuplicator.h"

#include "itkHistogramMatchingImageFilter.h"

#include "itkIdentityTransform.h"

#include "itkImage.h"
#include "itkImageDuplicator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

//#include "itkMultiplyByConstantImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkResampleImageFilter.h"

#include "itkTransformFactoryBase.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkVector.h"


// vnl
#include "vnl/vnl_matrix.h"


// local
#include "utilitiesImage.h"

namespace ShapeAnalysis
{
  /**
   * castItkImage
   */
  template< typename inputPixel_t, typename outputPixel_t > 
  typename itk::Image<outputPixel_t, 3 >::Pointer
  castItkImage( typename itk::Image<inputPixel_t, 3>::Pointer inputImage )
  {
    const unsigned int Dimension = 3;

    typedef itk::Image<inputPixel_t, Dimension> inputImage_t;
    typedef itk::Image<outputPixel_t, Dimension> outputImage_t;

    typedef itk::CastImageFilter< inputImage_t, outputImage_t > itkCastFilter_t;

    typename itkCastFilter_t::Pointer caster = itkCastFilter_t::New();
    caster->SetInput( inputImage );
    caster->Update();


    return caster->GetOutput();
  }


  template<typename image_t>
  double getVol(typename image_t::Pointer img, typename image_t::PixelType thld)
  {
    typedef itk::ImageRegionConstIterator<image_t> ImageRegionConstIterator_t;
    ImageRegionConstIterator_t it(img, img->GetLargestPossibleRegion() );

    double cell = (img->GetSpacing()[0])*(img->GetSpacing()[1])*(img->GetSpacing()[2]);

    double v = 0.0;

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        typename image_t::PixelType f = it.Get();

        v += (f>thld?cell:0.0);
      }

    return v;
  }



  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                 \
                  typename input_image_t::PixelType lowerT,                         \
                  typename input_image_t::PixelType upperT, \
                  typename output_image_t::PixelType insideValue,               \
                  typename output_image_t::PixelType outsideValue)
  {
    /**
     * O(x) :=    I(x) \in [lowerT, upperT] ? insideValue : outsideValue
     */

    //tst
    //   std::cout<<lowerT<<std::endl;
    //   std::cout<<upperT<<std::endl;
    //tst//

    typedef itk::BinaryThresholdImageFilter<input_image_t, output_image_t> binaryThresholdImageFilter_t;

    typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
    thlder->SetInput(input);
    thlder->SetInsideValue(insideValue);
    thlder->SetOutsideValue(outsideValue);
    thlder->SetUpperThreshold(upperT);
    thlder->SetLowerThreshold(lowerT);
    thlder->Update();
  
    return thlder->GetOutput();
  }


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                \
                  typename input_image_t::PixelType thld,               \
                  typename output_image_t::PixelType insideValue)
  {
    typename input_image_t::PixelType lowerT = thld;
    typename input_image_t::PixelType upperT = static_cast<typename input_image_t::PixelType>(1e16);
    typename output_image_t::PixelType outsideValue = 0;

    return binarilizeImage<input_image_t, output_image_t>(input, lowerT, upperT, insideValue, outsideValue);
  }


  template<typename ImageType>
  void binarilizeImageSeries(std::vector<typename ImageType::Pointer>& imgList)
  {
    std::size_t n = imgList.size();

    for (std::size_t it = 0; it < n; ++it)
      {
        typename ImageType::PixelType lowerThresholdInclusive = 1;
        typename ImageType::PixelType upperThresholdInclusive = std::numeric_limits<typename ImageType::PixelType>::max();

        imgList[it] = binarilizeImage<ImageType, ImageType>(imgList[it], \
                                                            lowerThresholdInclusive, \
                                                            upperThresholdInclusive, \
                                                            1,          \
                                                            0);
      }
  
    return;
  }


  template<typename image_t>
  typename image_t::RegionType
  computeNonZeroRegion(typename image_t::Pointer img)
  {
    /**
     * Given the img, compute the region where outside this region,
     * the image is all zero.
     *
     * The minx, y, z are initialized as sizeX, y, z; then, whenever
     * encounter an non-zero voxel, the minx, y, z are updated
     * accordingly. Similar for maxX, y, z except that they are
     * intialized to 0, 0, 0
     */
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;


    imageRegion_t entireRegion = img->GetLargestPossibleRegion();

    long minX = entireRegion.GetSize()[0];
    long minY = entireRegion.GetSize()[1];
    long minZ = entireRegion.GetSize()[2];

    long maxX = 0;
    long maxY = 0;
    long maxZ = 0;

    //    std::cout<<"hahaha = "<<minX<<'\t'<<minY<<'\t'<<minZ<<'\t'<<maxX<<'\t'<<maxX<<'\t'<<maxX<<'\n';

    typedef itk::ImageRegionConstIteratorWithIndex< image_t > itkImageRegionConstIteratorWithIndex_t;

    itkImageRegionConstIteratorWithIndex_t it(img, entireRegion);

    char foundNonZero = 0;

    {
      imageIndex_t idx;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
          if (it.Get() != 0)
            {
              foundNonZero = 1;

              idx = it.GetIndex();

              minX = minX<idx[0]?minX:idx[0];
              minY = minY<idx[1]?minY:idx[1];
              minZ = minZ<idx[2]?minZ:idx[2];

              maxX = maxX>idx[0]?maxX:idx[0];
              maxY = maxY>idx[1]?maxY:idx[1];
              maxZ = maxZ>idx[2]?maxZ:idx[2];
            }
        }
    }

    imageRegion_t nonZeroRegion;

    if (1 == foundNonZero)
      {
        imageIndex_t startIdx;
        startIdx[0] = minX;
        startIdx[1] = minY;
        startIdx[2] = minZ;

        imageSize_t size;
        size[0] = maxX - minX;
        size[1] = maxY - minY;
        size[2] = maxZ - minZ;

        nonZeroRegion.SetSize( size );
        nonZeroRegion.SetIndex( startIdx );
      }
    else
      {
        imageIndex_t startIdx;
        startIdx[0] = 0;
        startIdx[1] = 0;
        startIdx[2] = 0;

        imageSize_t size;
        size[0] = entireRegion.GetSize()[0];
        size[1] = entireRegion.GetSize()[1];
        size[2] = entireRegion.GetSize()[2];

        nonZeroRegion.SetSize( size );
        nonZeroRegion.SetIndex( startIdx );
      }

    
    return nonZeroRegion;
  }


  /**
   * Enlarge the region by 1/5 at each end, care is taken at the
   * boundary.
   */
  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegion(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion)
  {
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;

    imageRegion_t entireRegion = img->GetLargestPossibleRegion();
    imageSize_t entireSize = entireRegion.GetSize();

    imageIndex_t start = nonZeroRegion.GetIndex();
    imageSize_t sz = nonZeroRegion.GetSize();

    start[0] = std::max(0l, static_cast<long>(start[0] - sz[0]/5));
    start[1] = std::max(0l, static_cast<long>(start[1] - sz[1]/5));
    start[2] = std::max(0l, static_cast<long>(start[2] - sz[2]/5));

    sz[0] = std::min(entireSize[0] - start[0], 7*sz[0]/5);
    sz[1] = std::min(entireSize[1] - start[1], 7*sz[1]/5);
    sz[2] = std::min(entireSize[2] - start[2], 7*sz[2]/5);

    
    /**********************************************************************************
    {
      //tst
      std::cout<<"\t\t start =    "<<start<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize =    "<<entireSize<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize[1] - start[1], 7*sz[1]/5   "<<entireSize[1] - start[1]<<'\t'<<7*sz[1]/5<<'\t'<<sz[1]<<std::endl<<std::flush;
      //tst//
    }
    **********************************************************************************/

    imageRegion_t largerRegion;
    largerRegion.SetSize( sz );
    largerRegion.SetIndex( start );

    return largerRegion;
  }


  /**
   * Extract the ROI from the image using the region
   */
  template<typename image_t>
  typename image_t::Pointer
  extractROI(typename image_t::Pointer img, typename image_t::RegionType region)
  {
    typedef itk::RegionOfInterestImageFilter<image_t, image_t> itkRegionOfInterestImageFilter_t;

    typename itkRegionOfInterestImageFilter_t::Pointer ROIfilter = itkRegionOfInterestImageFilter_t::New();
    ROIfilter->SetInput( img );
    ROIfilter->SetRegionOfInterest( region );
    ROIfilter->Update();

    return ROIfilter->GetOutput();
  }


  /**
   * Crop the image by the non-zero region given by the mask 
   */
  /**
   * Crop the mask by its non-zero region
   */
  template<typename MaskImageType >
  typename MaskImageType::Pointer
  cropNonZeroRegionFromImage(typename MaskImageType::Pointer mask)
  {
    typename MaskImageType::RegionType ROIRegion = computeNonZeroRegion<MaskImageType>(mask);

    typename MaskImageType::RegionType enlargedROIRegion = enlargeNonZeroRegion<MaskImageType>(mask, ROIRegion);

    typename MaskImageType::Pointer ROIMask = extractROI<MaskImageType>(mask, enlargedROIRegion);

    return ROIMask;
  }




  /** 
   * Generate an all-zero image the same size/origin/spacing/etc. as
   * referenceImg, inside of whick, the roiRegion is the roiImg
   */
  template<typename image_t>
  typename image_t::Pointer
  antiExtractROI(typename image_t::ConstPointer roiImg, const typename image_t::RegionType roiRegion, \
                 typename image_t::ConstPointer referenceImg)
  {

    /********************************************************************************
    {
      //tst
      std::cout<<"\t in antiExtractROI\n"<<std::flush;
      std::cout<<"\t roiImg.GetLargestPossibleRegion() = "<<roiImg->GetLargestPossibleRegion()<<std::endl<<std::flush;
      std::cout<<"\t roiRegion = "<<roiRegion<<std::endl<<std::flush;
      std::cout<<"\t referenceImg.GetLargestPossibleRegion() = "<<referenceImg->GetLargestPossibleRegion()<<std::endl<<std::flush;
      //tst//
    }
    ********************************************************************************/

    typedef typename image_t::Pointer imagePointer_t;

    imagePointer_t largeImage = image_t::New();
    largeImage->SetRegions( referenceImg->GetLargestPossibleRegion() );
    largeImage->Allocate();

    largeImage->FillBuffer(0);
    largeImage->CopyInformation(referenceImg);


    typedef itk::ImageRegionIterator< image_t > itkImageRegionIterator_t;
    typedef itk::ImageRegionConstIterator< image_t > itkImageRegionConstIterator_t;

    {
      itkImageRegionConstIterator_t itROI(roiImg, roiImg->GetLargestPossibleRegion());
      itkImageRegionIterator_t itNew(largeImage, roiRegion);

      itROI.GoToBegin();
      itNew.GoToBegin();
      for (; !itROI.IsAtEnd(); ++itROI, ++itNew)
        {
          itNew.Set(itROI.Get());
        }
    }

    return largeImage;
  }



  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegionByOnePixel(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion)
  {
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;

    imageRegion_t entireRegion = img->GetLargestPossibleRegion();
    imageSize_t entireSize = entireRegion.GetSize();

    imageIndex_t start = nonZeroRegion.GetIndex();
    imageSize_t sz = nonZeroRegion.GetSize();


    start[0] = std::max(0l, static_cast<long>(start[0] - 1));
    start[1] = std::max(0l, static_cast<long>(start[1] - 1));
    start[2] = std::max(0l, static_cast<long>(start[2] - 1));

    sz[0] = std::min(entireSize[0] - start[0], sz[0] + 2);
    sz[1] = std::min(entireSize[1] - start[1], sz[1] + 2);
    sz[2] = std::min(entireSize[2] - start[2], sz[2] + 2);


    /**********************************************************************************    
    {
      //tst
      std::cout<<"\t\t start =    "<<start<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize =    "<<entireSize<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize[1] - start[1], 7*sz[1]/5   "<<entireSize[1] - start[1]<<'\t'<<7*sz[1]/5<<'\t'<<sz[1]<<std::endl<<std::flush;
      //tst//
    }
    ********************************************************************************/

    imageRegion_t largerRegion;
    largerRegion.SetSize( sz );
    largerRegion.SetIndex( start );

    return largerRegion;
  }


  /**
   * convert between probability image and logit image 
   */
  template<typename TInputImage, typename TOutputImage>
  typename TOutputImage::Pointer
  probabilityImageToLogitImage(typename TInputImage::Pointer inputImage)
  {
    typedef TInputImage InputImageType;
    typedef TOutputImage OutputImageType;

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typename TInputImage::RegionType region = inputImage->GetLargestPossibleRegion();

    typename OutputImageType::Pointer logitImage = OutputImageType::New();
    logitImage->SetRegions(region);
    logitImage->CopyInformation(inputImage);
    logitImage->Allocate();
    logitImage->FillBuffer(0.0);

    typedef itk::ImageRegionConstIterator<InputImageType> ImageRegionConstIterator;
    ImageRegionConstIterator citer(inputImage, region);

    typedef itk::ImageRegionIterator<OutputImageType> ImageRegionIterator;
    ImageRegionIterator iter(logitImage, region);

    citer.GoToBegin();
    iter.GoToBegin();

    for (; !citer.IsAtEnd(); ++citer, ++iter)
      {
        InputPixelType v = citer.Get();
        if (v > 1.0 || v < 0.0)
          {
            std::cerr<<"v = "<<v<<", Error: need 0.0 <= v <= 1.0\n";
            abort();
          }

        OutputPixelType vv = static_cast<OutputPixelType>(vcl_log(v + vnl_math::eps) - vcl_log(1.0 - v + vnl_math::eps));

        iter.Set(vv);
      }

    return logitImage;
  }

  template<typename TInputImage, typename TOutputImage>
  typename TOutputImage::Pointer
  logitImageToProbabilityImage(typename TInputImage::Pointer inputImage)
  {
    typedef TInputImage InputImageType;
    typedef TOutputImage OutputImageType;

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typename TInputImage::RegionType region = inputImage->GetLargestPossibleRegion();

    typename OutputImageType::Pointer probImage = OutputImageType::New();
    probImage->SetRegions(region);
    probImage->CopyInformation(inputImage);
    probImage->Allocate();
    probImage->FillBuffer(0.0);

    typedef itk::ImageRegionConstIterator<OutputImageType> ImageRegionConstIterator;
    ImageRegionConstIterator citer(inputImage, region);

    typedef itk::ImageRegionIterator<TInputImage> ImageRegionIterator;
    ImageRegionIterator iter(probImage, region);

    citer.GoToBegin();
    iter.GoToBegin();

    for (; !citer.IsAtEnd(); ++citer, ++iter)
      {
        InputPixelType v = citer.Get();

        OutputPixelType vv = static_cast<OutputPixelType>(1.0 / (1.0 + vcl_exp(-v) ));

        iter.Set(vv);
      }

    return probImage;
  }



  template< typename ImageType, typename OutputImageType >
  typename OutputImageType::Pointer
  isotropicizeImage(typename ImageType::Pointer img)
  {
    typedef typename ImageType::SpacingType SpacingType;
    SpacingType inputSpacing = img->GetSpacing();
    float smallestSpacing = inputSpacing[0] > inputSpacing[1]?inputSpacing[1]:inputSpacing[0];
    smallestSpacing = smallestSpacing > inputSpacing[2]?inputSpacing[2]:smallestSpacing;

    SpacingType outputSpacing;
    outputSpacing.Fill(smallestSpacing);

    typedef typename ImageType::SizeType SizeType;
    SizeType inputSize = img->GetLargestPossibleRegion().GetSize();
    SizeType outputSize;
    outputSize[0] = static_cast<typename SizeType::SizeValueType>(inputSize[0]*(inputSpacing[0]/outputSpacing[0]));
    outputSize[1] = static_cast<typename SizeType::SizeValueType>(inputSize[1]*(inputSpacing[1]/outputSpacing[1]));
    outputSize[2] = static_cast<typename SizeType::SizeValueType>(inputSize[2]*(inputSpacing[2]/outputSpacing[2]));

    typedef itk::IdentityTransform<double, ImageType::ImageDimension> IDTransformType; // has to be double, float causes compile error, not sure why
    typedef itk::ResampleImageFilter<ImageType, OutputImageType> ResampleImageFilterType;
    typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
    resample->SetInput(img);
    resample->SetOutputOrigin(img->GetOrigin());
    resample->SetOutputDirection(img->GetDirection());
    resample->SetSize(outputSize);
    resample->SetOutputSpacing(outputSpacing);
    resample->SetTransform(IDTransformType::New());
    resample->UpdateLargestPossibleRegion();
 
    return resample->GetOutput();
  }


}// namespace ShapeAnalysis

#endif
