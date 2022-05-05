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


#ifndef utilitiesVTK_h_
#define utilitiesVTK_h_

#include <vector>
#include <string>


// itk
//#include "itkImage.h"

// vtk
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace ShapeAnalysis
{
  /********************************************************************************
   * Attach a scalar field to the points on the polydata
   */
  template< typename ScalarType >
  void
  appendScalarFieldToPolydataPoints(vtkSmartPointer<vtkPolyData> pd, const std::vector<ScalarType>& scalarField, \
                                    const std::string& name);

  template< typename ScalarType >
  void
  appendScalarFieldToPolydataPoints(vtkSmartPointer<vtkPolyData> pd, const ScalarType* scalarField, long n, const std::string& name);

  // /*******************************************************************************
  //  * extract iso surface from gray scale image, Adopted from
  //  * Slicer/Applications/CLI/GrayscaleModelMaker.cxx
  //  */
  // template<typename TNull>
  // vtkSmartPointer<vtkPolyData>
  // isoSurf(vtkSmartPointer<vtkImageData> image, float thld);

  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  polydataToggleLPSAndRAS(vtkSmartPointer<vtkPolyData> pd);


  /* Given a list of polydata with the same   */
  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  polydataPointDataMin(std::vector< vtkSmartPointer<vtkPolyData> > pdlist);




}// namespace ShapeAnalysis


#include "utilitiesVTK.hxx"

#endif
