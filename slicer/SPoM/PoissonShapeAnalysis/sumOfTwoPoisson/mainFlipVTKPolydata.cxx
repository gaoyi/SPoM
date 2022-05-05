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
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>

#include <vtkProperty.h>


// local
#include "utilitiesIO.h"

int main(int argc, char** argv)
{
  if (argc < 3)
    {
      std::cerr<<"args: input.vtp output.vtp\n";
      exit(-1);
    }

  std::string inputName(argv[1]);
  std::string outputName(argv[2]);

  vtkSmartPointer<vtkPolyData> pd = ShapeAnalysis::readPolyDataFromVtpFile<char>(inputName.c_str());

  vtkSmartPointer<vtkTransform> translation = vtkSmartPointer<vtkTransform>::New();
  translation->Scale(-1.0, -1.0, 1.0);


  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputData(pd);
  transformFilter->SetTransform(translation);
  transformFilter->Update();

  ShapeAnalysis::writePolydataToXMLFile<char>(outputName.c_str(), transformFilter->GetOutput());


  return EXIT_SUCCESS;
}
