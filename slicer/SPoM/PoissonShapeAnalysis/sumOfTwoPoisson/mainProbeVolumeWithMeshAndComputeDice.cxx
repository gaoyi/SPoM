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


#include <vector>


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkShortArray.h>


#include "itkContinuousIndex.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"


// local
#include "utilitiesIO.h"

int main(int argc, char** argv)
{
  if (argc < 3)
    {
      std::cerr<<"args: input.vtp inputVolume.nrrd RAS/LPS(1/-1)\n";
      exit(-1);
    }

  std::string inputName(argv[1]);
  std::string inputVolumeName(argv[2]);

  int ras = 1;
  if (argc >= 4)
    {
      int v = atoi(argv[3]);

      if (v != 1 && v != -1)
        {
          std::cerr<<"Error: RAS = 1, LPS = -1. Other not supported.\n";
          abort();
        }

      ras = v;
    }

  char interpolationType = 0;

  vtkSmartPointer<vtkPolyData> polydata = ShapeAnalysis::readPolyDataFromVtpFile<char>(inputName.c_str());

  vtkSmartPointer<vtkDoubleArray> pValueArray = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("pValueFDR"));

  // Get the number of points in the polydata
  vtkIdType idNumPointsInFile = polydata->GetNumberOfPoints();

  std::cout<<idNumPointsInFile<<std::endl<<std::flush;

  std::vector<double> pValues(idNumPointsInFile);

  if(pValueArray)
    {
      for(int i = 0; i < idNumPointsInFile; i++)
        {
          //std::cout << "Got array." << std::endl;
          pValues[i] = pValueArray->GetValue(i);
          /*
          //if the array held arrays instead of scalars, you would use this:
          double location[3];
          array->GetTupleValue(i, location);
          std::cout << "Location: " << Location[0] << ","  << Location[1] << "," << Location[2] << std::endl;
          */
          //std::cout << "Distance: " << dist << std::endl;
          //std::cout << dist <<",";
        }
    }//end if(array)
  else
    {
      std::cout << "no pvalue array." << std::endl;
    }



  typedef short ImagePixelType;
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


  // vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
  // weights->SetNumberOfValues(polydata->GetNumberOfPoints());

  double numItersection = 0.0;
  double numUnion = 0.0;

  double inGroundTruth = 0.0;
  double inPValue = 0.0;

  itk::ContinuousIndex<double, 3> itkContinuousIndex;
  itk::Point< CoordinateRepresentationType, 3 > itkPoint;
  for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
      double p[3];
      polydata->GetPoint(i,p);
      // This is identical to:
      // polydata->GetPoints()->GetPoint(i,p);

      itkPoint[0] = ras*p[0];
      itkPoint[1] = ras*p[1];
      itkPoint[2] = p[2];

      image->TransformPhysicalPointToContinuousIndex(itkPoint, itkContinuousIndex);

      short l = 0.0;

      if (interpolator->IsInsideBuffer(itkContinuousIndex))
        {
          l = interpolator->EvaluateAtContinuousIndex(itkContinuousIndex);

          // if (l == 0)
          //   {
          //     itk::ContinuousIndex<double, 3> itkContinuousIndex1 = itkContinuousIndex;

          //     for (int itz = -1; itz <= 1; ++itz)
          //       {
          //         itkContinuousIndex1[2] = itkContinuousIndex[2] + itz;

          //         for (int ity = -1; ity <= 1; ++ity)
          //           {
          //             itkContinuousIndex1[1] = itkContinuousIndex[1] + ity;

          //             for (int itx = -1; itx <= 1; ++itx)
          //               {
          //                 itkContinuousIndex1[0] = itkContinuousIndex[0] + itx;

          //                 if (interpolator->IsInsideBuffer(itkContinuousIndex1))
          //                   {
          //                     short l1 = interpolator->EvaluateAtContinuousIndex(itkContinuousIndex1);

          //                     l = l1>l?l1:l;
          //                   }
          //               }
          //           }
          //       }
          //   }


          std::cout<<l<<"     "<<pValues[i]<<std::endl<<std::flush;

          if (2 == l)
            {
              ++inGroundTruth;
            }

          if (0.05 >= pValues[i])
            {
              ++inPValue;
            }


          if (2 == l || 0.05 >= pValues[i])
            {
              ++numUnion;
            }

          if (2 == l && 0.05 >= pValues[i])
            {
              ++numItersection;
            }
        }
    }

  std::cout<<numItersection/numUnion<<std::endl;
  std::cout<<(2.0*static_cast<double>(numItersection))/(static_cast<double>(inGroundTruth) + static_cast<double>(inPValue))<<std::endl;


  return EXIT_SUCCESS;
}
