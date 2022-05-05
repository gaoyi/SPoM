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


#ifndef ShapeSumOfTwoPoissonFilter_hxx_
#define ShapeSumOfTwoPoissonFilter_hxx_

#include <vector>

// itk
#include "itkImage.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

// local
#include "DirichletPoissonSolver3D.h"
#include "NeumannDirichletPoissonSolver3D.h"

#include "ShapeSumOfTwoPoissonFilter.h"


namespace ShapeAnalysis
{
  template< typename TShapeImageType >
  void
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::update()
  {
    //    _preprocess();

    _computeSumOfTwoPoisson();

    m_allDone = true;

    return;
  }

  template< typename TShapeImageType >
  void
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::_computeSumOfTwoPoisson()
  {
    /**
     * Compute Inside and outside Poissons. Order can NOT interchange!
     */

    //dbg
    std::cout<<"_computeInsidePoisson...\n"<<std::flush;
    //dbg, end
    _computeInsidePoisson();


    //dbg
    std::cout<<"_computeOutsidePoisson...\n"<<std::flush;
    //dbg, end
    _computeOutsidePoisson();

    return;
  }

  template< typename TShapeImageType >
  void
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::_computeInsidePoisson()
  {
    // set up Dirichelet BC mask. note that the boundary vovels are all 1
    m_DirichletBCMaskForPoissonSolver = CharImageType::New();
    m_DirichletBCMaskForPoissonSolver->SetRegions(m_inputShape->GetLargestPossibleRegion());
    m_DirichletBCMaskForPoissonSolver->Allocate();
    m_DirichletBCMaskForPoissonSolver->FillBuffer(1);
    m_DirichletBCMaskForPoissonSolver->CopyInformation(m_inputShape);
    {
      long nx = m_inputShape->GetLargestPossibleRegion().GetSize()[0];
      long ny = m_inputShape->GetLargestPossibleRegion().GetSize()[1];
      long nz = m_inputShape->GetLargestPossibleRegion().GetSize()[2];

      typename CharImageType::IndexType idx;
      for (long iz = 1; iz < nz-1; ++iz)
        {
          idx[2] = iz;
          for (long iy = 1; iy < ny-1; ++iy)
            {
              idx[1] = iy;
              for (long ix = 1; ix < nx-1; ++ix)
                {
                  idx[0] = ix;

                  // Because above I set all to 1, so here I need to
                  // "release" those non-BC positions.
                  if (m_inputShape->GetPixel(idx) != 0)
                    {
                      m_DirichletBCMaskForPoissonSolver->SetPixel(idx, 0);
                    }
                }
            }
        }
    }


    // set up Dirichelet BC values
    float edgeBC = 0.0;

    m_DirichletBCForPoissonSolver = FloatImageType::New();
    m_DirichletBCForPoissonSolver->SetRegions(m_inputShape->GetLargestPossibleRegion());
    m_DirichletBCForPoissonSolver->Allocate();
    m_DirichletBCForPoissonSolver->FillBuffer(edgeBC); // a large value for image boundary
    m_DirichletBCForPoissonSolver->CopyInformation(m_inputShape);

    /**
     * Now solve
     */
    typedef DirichletPoissonSolver3D<CharImageType, FloatImageType> DirichletPoissonSolver3DType;
    DirichletPoissonSolver3DType foo;

    foo.setBoundaryConditionMask(m_DirichletBCMaskForPoissonSolver);
    foo.setBoundaryConditionValue(m_DirichletBCForPoissonSolver);

    foo.update();

    m_poissonImage = foo.getPoissonImage();


    return;
  }


  template< typename TShapeImageType >
  void
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::_computeOutsidePoisson()
  {
    // set up Dirichelet BC mask. note that the boundary vovels are all 1
    m_DirichletBCMaskForPoissonSolver = CharImageType::New();
    m_DirichletBCMaskForPoissonSolver->SetRegions(m_inputShape->GetLargestPossibleRegion());
    m_DirichletBCMaskForPoissonSolver->Allocate();
    m_DirichletBCMaskForPoissonSolver->FillBuffer(0);
    m_DirichletBCMaskForPoissonSolver->CopyInformation(m_inputShape);
    {
      long nx = m_inputShape->GetLargestPossibleRegion().GetSize()[0];
      long ny = m_inputShape->GetLargestPossibleRegion().GetSize()[1];
      long nz = m_inputShape->GetLargestPossibleRegion().GetSize()[2];

      typename CharImageType::IndexType idx;
      for (long iz = 1; iz < nz-1; ++iz)
        {
          idx[2] = iz;
          for (long iy = 1; iy < ny-1; ++iy)
            {
              idx[1] = iy;
              for (long ix = 1; ix < nx-1; ++ix)
                {
                  idx[0] = ix;

                  if (m_inputShape->GetPixel(idx) != 0)
                    {
                      m_DirichletBCMaskForPoissonSolver->SetPixel(idx, 1);
                    }
                }
            }
        }
    }


    // set up Dirichelet BC values
    float edgeBC = 0.0;

    m_DirichletBCForPoissonSolver = FloatImageType::New();
    m_DirichletBCForPoissonSolver->SetRegions(m_inputShape->GetLargestPossibleRegion());
    m_DirichletBCForPoissonSolver->Allocate();
    m_DirichletBCForPoissonSolver->FillBuffer(edgeBC); // a large value for image boundary
    m_DirichletBCForPoissonSolver->CopyInformation(m_inputShape);


    /**
     * Now solve
     */
    typedef NeumannDirichletPoissonSolver3D<CharImageType, FloatImageType> NeumannDirichletPoissonSolver3DType;
    NeumannDirichletPoissonSolver3DType foo;

    foo.setBoundaryConditionMask(m_DirichletBCMaskForPoissonSolver);
    foo.setBoundaryConditionValue(m_DirichletBCForPoissonSolver);

    foo.update();

    typename NeumannDirichletPoissonSolver3DType::FloatImagePointer outsidePoisson = foo.getPoissonImage();


    writeImage<FloatImageType>(outsidePoisson, "outsidePoisson.nrrd");
    writeImage<FloatImageType>(m_poissonImage, "insidePoisson.nrrd");



    // now combine the inside and outside Poisson's
    {
      long nx = m_inputShape->GetLargestPossibleRegion().GetSize()[0];
      long ny = m_inputShape->GetLargestPossibleRegion().GetSize()[1];
      long nz = m_inputShape->GetLargestPossibleRegion().GetSize()[2];

      typename CharImageType::IndexType idx;
      for (long iz = 0; iz < nz; ++iz)
        {
          idx[2] = iz;
          for (long iy = 0; iy < ny; ++iy)
            {
              idx[1] = iy;
              for (long ix = 0; ix < nx; ++ix)
                {
                  idx[0] = ix;

                  if (m_inputShape->GetPixel(idx) == 0)
                    {
                      m_poissonImage->SetPixel(idx, -outsidePoisson->GetPixel(idx));
                    }
                }
            }
        }
    }


    return;
  }


  template< typename TShapeImageType >
  typename ShapeSumOfTwoPoissonFilter<TShapeImageType>::FloatImagePointer
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::getPoissonImage()
  {
    if (m_allDone)
      {
        return m_poissonImage;
      }
    else
      {
        std::cerr<<"Error: not done.\n";
        abort();
      }
  }

  template< typename TShapeImageType >
  void
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::setInputShape(ShapeImagePointer inputShape)
  {
    m_inputShape = inputShape;
    return;
  }

  // template< typename TShapeImageType >
  // void
  // ShapeSumOfTwoPoissonFilter<TShapeImageType>::setNumberOfIterations(long n)
  // {
  //   m_numIter = n;
  //   return;
  // }


  template< typename TShapeImageType >
  long
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::_numberOfNonZeroVoxels(ShapeImagePointer img)
  {
    typedef itk::ImageRegionConstIterator<ShapeImageType> ImageRegionConstIterator;
    ImageRegionConstIterator it(img, img->GetLargestPossibleRegion() );

    long v = 0;

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        typename ShapeImageType::PixelType f = it.Get();

        v += (f>0?1:0);
      }

    return v;
  }

    
  template< typename TShapeImageType >
  ShapeSumOfTwoPoissonFilter<TShapeImageType>::ShapeSumOfTwoPoissonFilter()
  {
    m_allDone = false;

    //m_numIter = 1000;

    m_dt = 0.1;
    
    return;
  }

}// ShapeAnalysis



#endif
