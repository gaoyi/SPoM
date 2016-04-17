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


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>


#include "itkContinuousIndex.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"


// local
#include "utilitiesIO.h"

int main(int argc, char** argv)
{
  if (argc < 4)
    {
      std::cerr<<"args: input.vtp inputVolume.nrrd output.vtp\n";
      exit(-1);
    }

  std::string inputName(argv[1]);
  std::string inputVolumeName(argv[2]);
  std::string outputName(argv[3]);

  char interpolationType = 0;

  vtkSmartPointer<vtkPolyData> polydata = ShapeAnalysis::readPolyDataFromVtpFile<char>(inputName.c_str());

  typedef double ImagePixelType;
  typedef itk::Image<ImagePixelType, 3> ImageType;

  ImageType::Pointer image = ShapeAnalysis::readImage<ImageType>(inputVolumeName.c_str());


  typedef double CoordinateRepresentationType ;

  itk::InterpolateImageFunction<ImageType, CoordinateRepresentationType>::Pointer interpolator;
  if (interpolationType == 0)
    {
      // NN interpolation
      typedef itk::NearestNeighborInterpolateImageFunction<ImageType, CoordinateRepresentationType> InterpolatorType;
      interpolator = InterpolatorType::New();
    }
  else if(interpolationType == 1)
    {
      // linear
      typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateRepresentationType> InterpolatorType;
      interpolator = InterpolatorType::New();
    }
  else if(interpolationType == 2)
    {
      // bspline
      typedef double CoefficientType;
      typedef itk::BSplineInterpolateImageFunction<ImageType, CoordinateRepresentationType, CoefficientType> InterpolatorType;
      interpolator = InterpolatorType::New();
    }
  interpolator->SetInputImage(image);


  vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
  weights->SetNumberOfValues(polydata->GetNumberOfPoints());

  itk::ContinuousIndex<double, 3> itkContinuousIndex;
  itk::Point< CoordinateRepresentationType, 3 > itkPoint;
  for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
      double p[3];
      polydata->GetPoint(i,p);
      // This is identical to:
      // polydata->GetPoints()->GetPoint(i,p);

      itkPoint[0] = -p[0];
      itkPoint[1] = -p[1];
      itkPoint[2] = p[2];

      image->TransformPhysicalPointToContinuousIndex(itkPoint, itkContinuousIndex);

      double v = 0.0;

      if (interpolator->IsInsideBuffer(itkContinuousIndex))
        {
          v = interpolator->EvaluateAtContinuousIndex(itkContinuousIndex);

          if (fabs(v) < 1e-5)
            {
              itk::ContinuousIndex<double, 3> itkContinuousIndex1 = itkContinuousIndex;

              for (int itz = -1; itz <= 1; ++itz)
                {
                  itkContinuousIndex1[2] = itkContinuousIndex[2] + itz;

                  for (int ity = -1; ity <= 1; ++ity)
                    {
                      itkContinuousIndex1[1] = itkContinuousIndex[1] + ity;

                      for (int itx = -1; itx <= 1; ++itx)
                        {
                          itkContinuousIndex1[0] = itkContinuousIndex[0] + itx;

                          if (interpolator->IsInsideBuffer(itkContinuousIndex1))
                            {
                              double v1 = interpolator->EvaluateAtContinuousIndex(itkContinuousIndex1);

                              v = v1>v?v1:v;
                            }
                        }
                    }
                }
            }
        }

      weights->SetValue(i, v);
      //std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
    }

  polydata->GetPointData()->SetScalars(weights);

  ShapeAnalysis::writePolydataToXMLFile<char>(outputName.c_str(), polydata);


  return EXIT_SUCCESS;
}
