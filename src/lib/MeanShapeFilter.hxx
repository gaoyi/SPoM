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


#ifndef MeanShapeFilter_hxx_
#define MeanShapeFilter_hxx_

#include <vector>
#include <cmath>

// itk
#include "itkAddImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"

#include "utilitiesImage.h"

#include "MeanShapeFilter.h"

// for dbg
#include "utilitiesIO.h"


namespace ShapeAnalysis
{
  template< typename TInputShapeType, typename TOutputShapeType >
  void
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::update()
  {
    if (m_meanShapeType == 0)
      {
        _computeLogoddsMeanShape();
      }
    else if (m_meanShapeType == 1)
      {
        _computeAverageBinaryMeanShape();
      }


    // //dbg
    // writeImage<FloatImageType>(m_floatMeanShape, "meanShapeFloat.nrrd");
    // //dbg, end


    // //dbg
    // writeImage<OutputShapeType>(m_meanShape, "meanShapeBin.nrrd");
    // //dbg, end


    m_allDone = true;

    return;
  }


  template< typename TInputShapeType, typename TOutputShapeType >
  void
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::setMeanShapeType(char meanShapeType)
  {
    if (meanShapeType > 1 || meanShapeType < 0)
      {
        std::cerr<<"Error. non-supported mean shape type.\n";
        std::cerr<<"       supported ones: 0. (default) log-odds mean shape; 1. average of binaries\n";
        abort();
      }
    else
      {
        m_meanShapeType = meanShapeType;
      }
    
    return;
  }

  template< typename TInputShapeType, typename TOutputShapeType >
  void
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::setThresholds(FloatType thldLower, FloatType thldUpper)
  {
    if (thldUpper < thldLower)
      {
        std::cerr<<"Error: thldUpper < thldLower\n";
        abort();
      }

    m_inputThresholdUpper = thldUpper;
    m_inputThresholdLower = thldLower;

    return;
  }

  template< typename TInputShapeType, typename TOutputShapeType >
  void
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::setInputShapes(const std::vector<InputShapePointer>& inputShapes)
  {
    RegionType region = inputShapes[0]->GetLargestPossibleRegion();

    std::size_t n = inputShapes.size();
    for (std::size_t it = 0; it < n; ++it)
      {
        if (inputShapes[it]->GetLargestPossibleRegion() != region)
          {
            std::cerr<<"Error: region of the "<<it<<"-th shape does not match the first shape.\n";
            abort();
          }
      }

    m_inputShapes = inputShapes;

    return;
  }


  template< typename TInputShapeType, typename TOutputShapeType >
  typename MeanShapeFilter<TInputShapeType, TOutputShapeType>::OutputShapePointer 
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::getMeanShape()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: mean shape not computed yet.\n";
        abort();
      }

    return m_meanShape;
  }

  template< typename TInputShapeType, typename TOutputShapeType >
  typename MeanShapeFilter<TInputShapeType, TOutputShapeType>::FloatImagePointer
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::getFloatMeanShape()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: mean shape not computed yet.\n";
        abort();
      }

    return m_floatMeanShape;
  }


  template< typename TInputShapeType, typename TOutputShapeType >
  void
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::_computeAverageBinaryMeanShape()
  {
    /**
     * 1. Get (floating typed) binary image from each input shape.
     * 2. average them
     * 3. threshold it
     *
     * I could have performed the above steps separately, however,
     * that will cause too many copies of shapes (floating typed) in
     * memory and i'm a little worried about that. So let me just do
     * everything in a single shot.
     */
    
    long n = static_cast<long>(m_inputShapes.size());
    FloatType inc = 1.0/static_cast<FloatType>(n);

    RegionType region = m_inputShapes[0]->GetLargestPossibleRegion();

    m_floatMeanShape = FloatImageType::New();
    m_floatMeanShape->SetRegions(region);
    m_floatMeanShape->Allocate();
    m_floatMeanShape->FillBuffer(0.0);
    m_floatMeanShape->CopyInformation(m_inputShapes[0]);

    typedef itk::ImageRegionIterator<FloatImageType> ImageRegionIterator;
    typedef itk::ImageRegionConstIterator<FloatImageType> ImageRegionConstIterator;

    for (long it = 0; it < n; ++it)
      {
        ImageRegionConstIterator citer(m_inputShapes[it], region);
        ImageRegionIterator iter(m_floatMeanShape, region);

        citer.GoToBegin();
        iter.GoToBegin();

        for (; !citer.IsAtEnd(); ++citer, ++iter)
          {
            FloatType v = static_cast<FloatType>(citer.Get());
            if (m_inputThresholdLower <= v && v <= m_inputThresholdUpper)
              {
                iter.Set(iter.Get() + inc);
              }
          }
      }


    m_meanShape = binarilizeImage<FloatImageType, OutputShapeType>(m_floatMeanShape, \
                                                                   0.5, itk::NumericTraits< FloatType >::max(), \
                                                                   static_cast<OutputPixelType>(1), \
                                                                   static_cast<OutputPixelType>(0));

    return;
  }

  template< typename TInputShapeType, typename TOutputShapeType >
  void
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::_computeLogoddsMeanShape()
  {
    /**
     * 1. Get (floating typed) binary image from each input shape.
     * 2. construct SDF
     * 3. logodds of SDf
     *
     * The itk SDF filter should be multi-threaded so no need to
     * parallel here. this is not the main time sink.
     */
    
    long n = static_cast<long>(m_inputShapes.size());
    FloatType inc = 1.0/static_cast<FloatType>(n);

    RegionType region = m_inputShapes[0]->GetLargestPossibleRegion();

    m_floatMeanShape = FloatImageType::New();
    m_floatMeanShape->SetRegions(region);
    m_floatMeanShape->Allocate();
    m_floatMeanShape->FillBuffer(0.0);
    m_floatMeanShape->CopyInformation(m_inputShapes[0]);

    typedef itk::ImageRegionIterator<FloatImageType> ImageRegionIterator;
    typedef itk::ImageRegionConstIterator<FloatImageType> ImageRegionConstIterator;

    typedef itk::SignedDanielssonDistanceMapImageFilter< FloatImageType, FloatImageType > SDFFilterType;

#pragma omp parallel for
    for (long it = 0; it < n; ++it)
      {
        FloatImagePointer binInput                                      \
          = binarilizeImage<InputShapeType, FloatImageType>(m_inputShapes[it], \
                                                            static_cast<InputPixelType>(m_inputThresholdLower), \
                                                            itk::NumericTraits< InputPixelType >::max(), \
                                                            1.0, 0.0);

        typename SDFFilterType::Pointer SDFFilter = SDFFilterType::New();
        SDFFilter->SetInput(binInput);
        SDFFilter->Update();

        FloatImagePointer sdf = SDFFilter->GetOutput();

#pragma omp critical
        {
          ImageRegionConstIterator citer(sdf, region);
          ImageRegionIterator iter(m_floatMeanShape, region);

          citer.GoToBegin();
          iter.GoToBegin();

          for (; !citer.IsAtEnd(); ++citer, ++iter)
            {
              FloatType phi = static_cast<FloatType>(citer.Get());

              FloatType logodds = 0.5*(1.0 + tanh(phi/m_epsilonForHeviside))*inc;

              iter.Set(iter.Get() + logodds);
            }
        }
      }


    m_meanShape = binarilizeImage<FloatImageType, OutputShapeType>(m_floatMeanShape, \
                                                                   0.5, itk::NumericTraits< FloatType >::max(), \
                                                                   static_cast<OutputPixelType>(0), \
                                                                   static_cast<OutputPixelType>(1));

    return;
  }



  template< typename TInputShapeType, typename TOutputShapeType >
  MeanShapeFilter<TInputShapeType, TOutputShapeType>::MeanShapeFilter()
  {
    // to inclose 1 for binary case
    m_inputThresholdUpper = itk::NumericTraits< FloatType >::max();
    m_inputThresholdLower = 0.5;

    m_meanShapeType = 0;

    m_allDone = false;

    m_epsilonForHeviside = 0.1;
    
    return;
  }


}// ShapeAnalysis



#endif
