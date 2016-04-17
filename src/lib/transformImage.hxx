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


#ifndef transformImage_hxx_
#define transformImage_hxx_


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// itk
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include "itkVersorRigid3DTransform.h"
#include "itkSimilarity3DTransform.h"

#include "itkWarpImageFilter.h"

// local
#include "transformImage.h"

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
                 typename OutputImageType::PixelType fillInValue, char interpolationType)
  {
    typedef double CoordinateRepresentationType ;

    typename itk::InterpolateImageFunction<InputImageType, CoordinateRepresentationType>::Pointer interpolator;

    if (interpolationType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, CoordinateRepresentationType> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<InputImageType, CoordinateRepresentationType> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 2)
      {
        // bspline
        typedef double CoefficientType;
        typedef itk::BSplineInterpolateImageFunction<InputImageType, CoordinateRepresentationType, CoefficientType> InterpolatorType;
        interpolator = InterpolatorType::New();
      }


    typedef itk::ResampleImageFilter< InputImageType, OutputImageType > ResampleFilterType;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    resampler->SetTransform( transform );
    resampler->SetInput( inputImage );
    resampler->SetSize( referenceImage->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin(  referenceImage->GetOrigin() );
    resampler->SetOutputSpacing( referenceImage->GetSpacing() );
    resampler->SetOutputDirection( referenceImage->GetDirection() );
    resampler->SetDefaultPixelValue( fillInValue );
    resampler->SetInterpolator(  interpolator  );
    resampler->Update();
  
    return resampler->GetOutput();
  }



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
                 typename OutputImageType::PixelType fillInValue, char interpolationType)
{
    typedef double CoordinateRepresentationType ;

    typename itk::InterpolateImageFunction<InputImageType, CoordinateRepresentationType>::Pointer interpolator;

    if (interpolationType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, CoordinateRepresentationType> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<InputImageType, CoordinateRepresentationType> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 2)
      {
        // bspline
        typedef double CoefficientType;
        typedef itk::BSplineInterpolateImageFunction<InputImageType, CoordinateRepresentationType, CoefficientType> InterpolatorType;
        interpolator = InterpolatorType::New();
      }


    typedef itk::ResampleImageFilter< InputImageType, OutputImageType > ResampleFilterType;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    resampler->SetTransform( transform );
    resampler->SetInput( inputImage );
    resampler->SetSize( referenceImage->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin(  referenceImage->GetOrigin() );
    resampler->SetOutputSpacing( referenceImage->GetSpacing() );
    resampler->SetOutputDirection( referenceImage->GetDirection() );
    resampler->SetDefaultPixelValue( fillInValue );
    resampler->SetInterpolator(  interpolator  );
    resampler->Update();
  
    return resampler->GetOutput();
  }



  /**
   * Warp image using vector image
   * 
   * interpolationType = 0 for NN interp, 1 for linear interp, 2 for
   * bspline interp. Default is 1
   */
  template<typename ReferenceImageType, typename MovingImageType, typename DisplacementFieldType>
  typename ReferenceImageType::Pointer
  warpImage(typename ReferenceImageType::Pointer img,                     \
            typename MovingImageType::Pointer movingImg,                  \
            typename DisplacementFieldType::Pointer displacementField, \
            typename ReferenceImageType::PixelType fillInVal, \
            char interpolationType)
  {
    typename itk::InterpolateImageFunction<MovingImageType, double>::Pointer interpolator;

    if (interpolationType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<MovingImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<MovingImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 2)
      {
        // bspline
        typedef itk::BSplineInterpolateImageFunction<MovingImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }


    typedef itk::WarpImageFilter<MovingImageType, ReferenceImageType, DisplacementFieldType > WarperType;

    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput( movingImg );
    warper->SetInterpolator( interpolator );
    warper->SetOutputSpacing( img->GetSpacing() );
    warper->SetOutputOrigin( img->GetOrigin() );
    warper->SetOutputDirection( img->GetDirection() );
    warper->SetEdgePaddingValue( fillInVal );
    warper->SetDisplacementField( displacementField );
    warper->Update();

    return warper->GetOutput();
  }





  template<typename InputImageType, typename OutputImageType, typename DeformationFieldType>
  typename OutputImageType::Pointer
  warpImage(typename DeformationFieldType::Pointer vectorField, \
            typename InputImageType::Pointer inputImage, typename OutputImageType::Pointer referenceImage, \
            typename OutputImageType::PixelType fillInValue, char interpolationType)
  {
    typename itk::InterpolateImageFunction<InputImageType, double>::Pointer interpolator;

    if (interpolationType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<InputImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpolationType == 2)
      {
        // bspline
        typedef itk::BSplineInterpolateImageFunction<InputImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }

    typedef itk::WarpImageFilter<InputImageType, OutputImageType, DeformationFieldType > WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput( inputImage );
    warper->SetInterpolator( interpolator );
    warper->SetOutputSpacing( referenceImage->GetSpacing() );
    warper->SetOutputOrigin( referenceImage->GetOrigin() );
    warper->SetDeformationField( vectorField );
    warper->Update();

    //tst
    //writeVectorImage<DeformationFieldType>(filter->GetOutput(), "defomationField.nrrd", 2);
    //tst//

    // InternalImageType::Pointer movingImgInternal = castItkImage<typename moving_image_t::PixelType, typename OutputImageType::PixelType>(movingImg);


    // typedef itk::CastImageFilter< InternalImageType, output_image_t > CastFilterType;
    // typename CastFilterType::Pointer  caster =  CastFilterType::New();
    // caster->SetInput( warper->GetOutput() );
    // caster->Update();

    return warper->GetOutput();

  }


}// namespce ShapeAnalysis



#endif //transformImage_hxx_
