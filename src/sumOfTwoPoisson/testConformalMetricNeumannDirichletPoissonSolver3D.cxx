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


#include <iostream>

#include "utilitiesIO.h"

#include "ConformalMetricNeumannDirichletPoissonSolver3D.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


template<typename LabelImageType, typename BoundaryConditionValueImageType>
typename BoundaryConditionValueImageType::Pointer
computeBoundaryConditionValueImage(const LabelImageType* labelImage);

int main(int argc, char** argv)
{
  if (argc < 5)
    {
      std::cerr<<"args: inputImage seedLabelImage outputImage beta\n";
      exit(-1);
    }

  std::string inputImageName(argv[1]);
  std::string seedLabelImageName(argv[2]);
  std::string outputImageName(argv[3]);
  double beta = std::atof(argv[4]);

  const unsigned int ImageDim = 3;
  typedef short InputImagePixelType;
  typedef itk::Image<InputImagePixelType, ImageDim> InputImageType;

  typedef short DirchletBoundaryConditionMaskPixelType;
  typedef itk::Image<DirchletBoundaryConditionMaskPixelType, ImageDim> DirchletBoundaryConditionMaskImageType;

  typedef double BoundaryValuePixelType;
  typedef itk::Image<BoundaryValuePixelType, ImageDim> BoundaryValueImageType;

  // Read in training images' names
  InputImageType::Pointer inputImage = ShapeAnalysis::readImage<InputImageType>(inputImageName.c_str());
  DirchletBoundaryConditionMaskImageType::Pointer seedLabelImage = ShapeAnalysis::readImage<DirchletBoundaryConditionMaskImageType>(seedLabelImageName.c_str());
  BoundaryValueImageType::Pointer bcvImage = computeBoundaryConditionValueImage<DirchletBoundaryConditionMaskImageType, BoundaryValueImageType>(seedLabelImage);

  typedef ShapeAnalysis::ConformalMetricNeumannDirichletPoissonSolver3D<InputImageType, DirchletBoundaryConditionMaskImageType, BoundaryValueImageType> ConformalMetricNeumannDirichletPoissonSolver3DType;
  ConformalMetricNeumannDirichletPoissonSolver3DType foo;

  foo.setInputImage(inputImage);
  foo.setBoundaryConditionMask(seedLabelImage);
  foo.setBoundaryConditionValue(bcvImage);
  foo.setBeta(beta);

  foo.update();

  ShapeAnalysis::writeImage<ConformalMetricNeumannDirichletPoissonSolver3DType::FloatImageType>(foo.getPoissonImage(), outputImageName.c_str());

  return 0;
}


template<typename LabelImageType, typename BoundaryConditionValueImageType>
typename BoundaryConditionValueImageType::Pointer
computeBoundaryConditionValueImage(const LabelImageType* labelImage)
{
  typename BoundaryConditionValueImageType::Pointer bcv = BoundaryConditionValueImageType::New();
  bcv->SetRegions(labelImage->GetLargestPossibleRegion());
  bcv->Allocate();
  bcv->CopyInformation(labelImage);
  bcv->FillBuffer(0.0);

  typedef itk::ImageRegionConstIterator<LabelImageType> ImageRegionConstIterator;
  ImageRegionConstIterator bciter(labelImage, labelImage->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<BoundaryConditionValueImageType> ImageRegionIterator;
  ImageRegionIterator iter(bcv, bcv->GetLargestPossibleRegion());

  bciter.GoToBegin();
  iter.GoToBegin();
  for (; !iter.IsAtEnd(); ++iter, ++bciter)
    {
      if (bciter.Get() == 1)
        {
          iter.Set(1.0);
        }
    }

  return bcv;
}
