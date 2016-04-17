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


/**
 * This defines the concreate itk mesh type which can be included in other applications.
 */

#ifndef meshTypedef_h_
#define meshTypedef_h_

#include "itkMesh.h"

namespace ShapeAnalysis
{
  ////////////////////////////////////////////////////////////////////////////////
  // Definitions
  // This define is needed to deal with double/float changes in VTK
#ifndef vtkFloatingPointType
#define vtkFloatingPointType double
#endif

  const unsigned int pointDimension   = 3;
  const unsigned int maxCellDimension = 3;

  typedef itk::Point<vtkFloatingPointType, pointDimension> ItkPoint;
  
  typedef itk::DefaultStaticMeshTraits< vtkFloatingPointType, pointDimension, maxCellDimension, vtkFloatingPointType, vtkFloatingPointType > MeshTraits;

  typedef itk::Mesh<vtkFloatingPointType, pointDimension, MeshTraits > MeshType;
}

#endif // meshTypedef_h_
