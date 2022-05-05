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

// system
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
// ITK Headers
#include "itkPoint.h"
#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkTriangleCell.h"

////////////////////////////////////////////////////////////////////////////////
// VTK headers
#include "vtkPoints.h"
#include "vtkCell.h"

#include "vtkPointData.h"
#include "vtkCellData.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include "vtkPolyDataNormals.h"

#include "vtkCellArray.h"

#include "vtkDataArray.h"
#include "vtkDataSetAttributes.h"

#include "vtkDoubleArray.h"

// local
#include "itkVtkMeshConversion.h"


namespace ShapeAnalysis
{
  MeshType::Pointer vtkPolyDataToITKMesh( vtkSmartPointer<vtkPolyData> polyData)
  {
    /// Create a new mesh
    MeshType::Pointer mesh = MeshType::New();

    /// Transfer the points from the vtkPolyData into the itk::Mesh
    const unsigned int numberOfPoints = polyData->GetNumberOfPoints();

    vtkSmartPointer<vtkPoints> vtkpoints = polyData->GetPoints();
  
    mesh->GetPoints()->Reserve( numberOfPoints );
  
    for(unsigned int p =0; p < numberOfPoints; p++)
      {
        double * apoint = vtkpoints->GetPoint( p );
        mesh->SetPoint( p, MeshType::PointType( apoint ));
      }

  
    /// Transfer the cells from the vtkPolyData into the itk::Mesh
    vtkSmartPointer<vtkCellArray> triangleStrips = polyData->GetStrips();


    const vtkIdType* cellPoints;
    vtkIdType numberOfCellPoints;

    /// First count the total number of triangles from all the triangle strips.
    unsigned int numberOfTriangles = 0;

    triangleStrips->InitTraversal();

    while( triangleStrips->GetNextCell( numberOfCellPoints, cellPoints ) )
      {
        numberOfTriangles += numberOfCellPoints-2;
      }

    vtkSmartPointer<vtkCellArray> polygons = polyData->GetPolys();

    polygons->InitTraversal();

    while( polygons->GetNextCell( numberOfCellPoints, cellPoints ) )
      {
        if( numberOfCellPoints == 3 )
          {
            numberOfTriangles ++;
          }
      }


    /// Reserve memory in the itk::Mesh for all those triangles
    mesh->GetCells()->Reserve( numberOfTriangles );

  
    /// Copy the triangles from vtkPolyData into the itk::Mesh
    typedef MeshType::CellType CellType;
    typedef itk::TriangleCell< CellType > TriangleCellType;

    int cellId = 0;

    /// first copy the triangle strips
    triangleStrips->InitTraversal();
    while( triangleStrips->GetNextCell( numberOfCellPoints, cellPoints ) )
      {
    
        unsigned int numberOfTrianglesInStrip = numberOfCellPoints - 2;

        unsigned long pointIds[3];
        pointIds[0] = cellPoints[0];
        pointIds[1] = cellPoints[1];
        pointIds[2] = cellPoints[2];

        for( unsigned int t=0; t < numberOfTrianglesInStrip; t++ )
          {
            MeshType::CellAutoPointer c;
            TriangleCellType * tcell = new TriangleCellType;
            tcell->SetPointIds( pointIds );
            c.TakeOwnership( tcell );
            mesh->SetCell( cellId, c );
            cellId++;
            pointIds[0] = pointIds[1];
            pointIds[1] = pointIds[2];
            pointIds[2] = cellPoints[t+3];
          }
      }

    /// then copy the normal triangles
    polygons->InitTraversal();
    while( polygons->GetNextCell( numberOfCellPoints, cellPoints ) )
      {
        if( numberOfCellPoints !=3 ) // skip any non-triangle.
          {
            continue;
          }
        MeshType::CellAutoPointer c;
        TriangleCellType * t = new TriangleCellType;
        //t->SetPointIds( static_cast<unsigned long*>(cellPoints) );
        t->SetPointIds( (unsigned long*)cellPoints );
        c.TakeOwnership( t );
        mesh->SetCell( cellId, c );
        cellId++;
      }

    return mesh;
  }

  vtkSmartPointer<vtkPolyData> ITKMeshToVtkPolyData( MeshType::Pointer mesh)
  {
    /// Creat a new vtkPolyData
    vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();

    /// Creat vtkPoints for insertion into newPolyData
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    /// Copy all points into the vtkPolyData structure
    typedef MeshType::PointsContainer::ConstIterator PointConstIterator;
    PointConstIterator pntIterator = mesh->GetPoints()->Begin();
    PointConstIterator pntItEnd = mesh->GetPoints()->End();

    for (int i = 0; pntIterator != pntItEnd; ++i, ++pntIterator)
      {
        ItkPoint pnt = pntIterator.Value();
        points->InsertPoint(i, pnt[0], pnt[1], pnt[2]);

        // //tst
        // std::cout<<i<<"-th point:  ";
        // std::cout<<pnt[0]<<std::endl;
        // std::cout<<"               "<<pntIterator.Value()<<std::endl;
        // //tst//
      }

    newPolyData->SetPoints(points);

    /// Copy all cells into the vtkPolyData structure
    /// Creat vtkCellArray into which the cells are copied
    vtkSmartPointer<vtkCellArray> triangle = vtkSmartPointer<vtkCellArray>::New();
  
    typedef MeshType::CellsContainer::ConstIterator CellConstIterator;
    CellConstIterator cellIt = mesh->GetCells()->Begin();
    CellConstIterator cellItEnd = mesh->GetCells()->End();


    {
      typedef MeshType::CellType CellType;

      for (int it = 0; cellIt != cellItEnd; ++it, ++cellIt)
        {
          CellType* cellptr = cellIt.Value();
          //    LineType * line = dynamic_cast<LineType *>( cellptr );
          //    std::cout << line->GetNumberOfPoints() << std::endl;
          //      std::cout << cellptr->GetNumberOfPoints() << std::endl;

          CellType::PointIdIterator pntIdIter = cellptr->PointIdsBegin();
          CellType::PointIdIterator pntIdEnd = cellptr->PointIdsEnd();
          vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

          for (; pntIdIter != pntIdEnd; ++pntIdIter)
            {
              pts->InsertNextId( *pntIdIter );

              // //tst
              // std::cout<<"           "<<pts[1]<<std::endl;
              // //tst//
            }
          triangle->InsertNextCell(pts);
        }
    }

    newPolyData->SetPolys(triangle);

    // //tst
    // std::cout<<"----------------------------------------here----------------------------------------\n"<<std::flush;
    // //tst//

    return newPolyData;
  }

} //namespace
////////////////////////////////////////////////////////////
