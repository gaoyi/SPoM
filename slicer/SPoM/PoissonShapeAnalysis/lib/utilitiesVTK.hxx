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


#ifndef utilitiesVTK_hxx_
#define utilitiesVTK_hxx_

#include <vector>
#include <string>

// vtk
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

// local
#include "utilitiesVTK.h"

namespace ShapeAnalysis
{
  /********************************************************************************
   * Attach a scalar field to the points on the polydata
   */
  template< typename ScalarType >
  void
  appendScalarFieldToPolydataPoints(vtkSmartPointer<vtkPolyData> pd, const std::vector<ScalarType>& scalarField, \
                                    const std::string& name)
  {
    appendScalarFieldToPolydataPoints<ScalarType>(pd, &scalarField[0], scalarField.size(), name);
      
    return;
  }


  template< typename ScalarType >
  void
  appendScalarFieldToPolydataPoints(vtkSmartPointer<vtkPolyData> pd, const ScalarType* scalarField, long n, const std::string& name)
  {
    if (n != pd->GetNumberOfPoints())
      {
        std::cerr<<"n != pd->GetNumberOfPoints()\n"<<std::endl;
        abort();
      }

    vtkSmartPointer<vtkDoubleArray> s = vtkSmartPointer<vtkDoubleArray>::New();
    s->SetName(name.c_str());
    s->SetNumberOfValues(n);
    for(vtkIdType i = 0; i < n; ++i)
      {
        s->SetValue(i, static_cast<double>(scalarField[i]));
      }

    //pd->GetPointData()->SetScalars(s);
    pd->GetPointData()->AddArray(s);

    return;
  }

  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  polydataToggleLPSAndRAS(vtkSmartPointer<vtkPolyData> oldpd)
  {
      vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
      pd->DeepCopy(oldpd);

    vtkIdType np = pd->GetNumberOfPoints();
    double p[3];

    for (vtkIdType ip = 0; ip < np; ++ip)
    {

        pd->GetPoint(ip, p);
        pd->GetPoints()->SetPoint(ip, -p[0], -p[1], p[2]);
        // This is identical to:
        // polydata->GetPoints()->GetPoint(i,p);
//        std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
    }

    return pd;
  }


  // /*******************************************************************************
  //  * extract iso surface from gray scale image, Adopted from
  //  * Slicer/Applications/CLI/GrayscaleModelMaker.cxx
  //  */
  // template<typename TNull>
  // vtkSmartPointer<vtkPolyData>
  // isoSurf(vtkSmartPointer<vtkImageData> image, float thld)
  // {
  //   vtkSmartPointer<vtkImageChangeInformation> ici = vtkSmartPointer<vtkImageChangeInformation>::New();
  //   ici->SetInput (image);
  //   ici->SetOutputSpacing( 1, 1, 1 );
  //   ici->SetOutputOrigin( 0, 0, 0 );
  //   ici->Update();

  //   vtkSmartPointer<vtkImageData> newImage = ici->GetOutput();
  //   newImage->Update();

  //   // Get the dimensions, marching cubes only works on 3d
  //   int extents[6];
  //   newImage->GetExtent(extents);
  //   if (debug)
  //     {
  //       std::cout << "Image data extents: " << extents[0] << " " << extents[1] << " " << extents[2] << " " << extents[3] << " " << extents[4] << " " << extents[5] << endl;
  //     }
  //   if (extents[0] == extents[1] ||
  //       extents[2] == extents[3] ||
  //       extents[4] == extents[5]) 
  //     {
  //       std::cerr << "The volume is not 3D.\n";
  //       std::cerr << "\tImage data extents: " << extents[0] << " " << extents[1] << " " << extents[2] << " " << extents[3] << " " << extents[4] << " " << extents[5] << endl;
  //       return EXIT_FAILURE;
  //     }
  //   // Get the RAS to IJK matrix and invert it to get the IJK to RAS which will need
  //   // to be applied to the model as it will be built in pixel space

  //   vtkSmartPointer<vtkTransform> transformIJKtoRAS = vtkSmartPointer<vtkTransform>::New();
  //   transformIJKtoRAS->SetMatrix(reader->GetRasToIjkMatrix());
  //   if (debug)
  //     {
  //       std::cout << "RasToIjk matrix from file = ";
  //       transformIJKtoRAS->GetMatrix()->Print(std::cout);
  //     }
  //   transformIJKtoRAS->Inverse();
  //   mcubes = vtkMarchingCubes::New();

  //   mcubes->SetInput(ici->GetOutput());
  //   mcubes->SetValue(0,Threshold);
  //   mcubes->ComputeScalarsOff();
  //   mcubes->ComputeGradientsOff();
  //   mcubes->ComputeNormalsOff();
  //   (mcubes->GetOutput())->ReleaseDataFlagOn();
  //   mcubes->Update();

  //   if (debug)
  //     {
  //       std::cout << "Number of polygons = " << (mcubes->GetOutput())->GetNumberOfPolys() << endl;
  //     }
    
  //   if ((transformIJKtoRAS->GetMatrix())->Determinant() < 0) 
  //     {
  //       if (debug)
  //         {
  //           std::cout << "Determinant " << (transformIJKtoRAS->GetMatrix())->Determinant() << " is less than zero, reversing...\n";
  //         }
  //       reverser = vtkReverseSense::New();
  //       reverser->SetInput(mcubes->GetOutput());
  //       reverser->ReverseNormalsOn();
  //       (reverser->GetOutput())->ReleaseDataFlagOn();
  //       // TODO: add progress
  //     }

  //   transformer = vtkTransformPolyDataFilter::New();
  //   if ((transformIJKtoRAS->GetMatrix())->Determinant() < 0) 
  //     {
  //       transformer->SetInput(reverser->GetOutput());
  //     } 
  //   else 
  //     {
  //       transformer->SetInput(mcubes->GetOutput());
  //     }

  //   transformer->SetTransform(transformIJKtoRAS);
  //   if (debug)
  //     {
  //       std::cout << "Transforming using inversed matrix:\n";
  //       transformIJKtoRAS->GetMatrix()->Print(std::cout);
  //     }

  //   // TODO: add progress
  //   (transformer->GetOutput())->ReleaseDataFlagOn();
    
    
  //   writer = vtkPolyDataWriter::New();
  //   writer->SetInput(transformer->GetOutput());
  //   writer->SetFileName(OutputGeometry.c_str());
  //   // TODO: add progress
  //   writer->Write();

  //   // Cleanup
  //   if (reader)
  //     {
  //       reader->Delete();
  //     }
  //   if (ici)
  //     {
  //       ici->Delete();
  //     }
  //   if (transformIJKtoRAS)
  //     {
  //       transformIJKtoRAS->Delete();
  //     }
  //   if (mcubes)
  //     {
  //       mcubes->Delete();
  //     }
  //   if (reverser)
  //     {
  //       reverser->Delete();
  //     }
  //   if (transformer)
  //     {
  //       transformer->Delete();
  //     }
  //   if (writer)
  //     {
  //       writer->Delete();
  //     }
  //   return EXIT_SUCCESS;

  // }
  

}// ShapeAnalysis



#endif
