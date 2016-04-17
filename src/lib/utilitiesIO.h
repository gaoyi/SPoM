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


#ifndef utilitiesIO_h_
#define utilitiesIO_h_


#include <string>

// itk
#include "itkAffineTransform.h"
#include "itkImage.h"

// vnl
#include "vnl/vnl_matrix.h"

// vtk
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

namespace ShapeAnalysis
{
  /**********************************************************************************
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName);

  /************************************************************************************
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName);

  /************************************************************************************
   * Read a series of images.
   */
  template< typename itkImage_t > 
  std::vector< typename itkImage_t::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList );

  /************************************************************************************
   * readTextLineToListOfString   
   */
  template<typename TNull>
  std::vector< std::string > readTextLineToListOfString(const char* textFileName);


  /************************************************************************************
   * write a component of a vector image
   */
  template< typename itkVectorImage_t > 
  void 
  writeVectorImage(typename itkVectorImage_t::Pointer img, const char *fileName, int component);


  /********************************************************************************
   * write vtkpolydata 
   */
  template< typename TNull >
  void writePolydataToXMLFile(const char *fileName, vtkSmartPointer<vtkPolyData> pd);

  /********************************************************************************
   * read vtkpolydata, from legacy .vtk file
   */
  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  readPolyDataFromVtkFile(const char *fileName);


  /********************************************************************************
   * read vtkpolydata, from .vtp file
   */
  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  readPolyDataFromVtpFile(const char *fileName);


  template< typename TNull >
  std::vector< vtkSmartPointer<vtkPolyData> >
  readPolyDataSeries(const std::vector<std::string>& names);
  

}// ShapeAnalysis


#include "utilitiesIO.hxx"

#endif
