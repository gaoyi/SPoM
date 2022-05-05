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


// local
#include "utilitiesIO.h"

int main(int argc, char** argv)
{
  if (argc < 4)
    {
      std::cerr<<"args: input.vtp value1 value2\n";
      exit(-1);
    }

  // In input.vtp, there are two scalar fields, f1 and f2. Then, for
  // all the vertex with f1==label1, say there are M of them. Then,
  // find how many (N) of them have f2==label2. Then output N/M

  std::string inputName(argv[1]);
  float v1 = atof(argv[2]);
  float v2 = atof(argv[3]);

  vtkSmartPointer<vtkPolyData> polydata = ShapeAnalysis::readPolyDataFromVtpFile<char>(inputName.c_str());

  // std::cout<<polydata->GetPointData()->GetArray("Scalars_")<<std::endl;
  // std::cout<<polydata->GetPointData()->GetArray("haha")<<std::endl;
  // std::cout<<polydata->GetPointData()->GetArray("pValueFDR")<<std::endl;

  vtkSmartPointer<vtkDoubleArray> labelArray = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("Scalars_"));

  // std::cout<<labelArray->GetNumberOfComponents()<<std::endl;
  // std::cout<<labelArray->GetNumberOfTuples()<<std::endl;

  vtkSmartPointer<vtkDoubleArray> pvalueArray = vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("pValueFDR"));
  // std::cout<<pvalueArray->GetNumberOfComponents()<<std::endl;
  // std::cout<<pvalueArray->GetNumberOfTuples()<<std::endl;


  if (labelArray && pvalueArray)
    {
      int n = labelArray->GetNumberOfTuples();
      if ( pvalueArray->GetNumberOfTuples() != n)
        {
          std::cerr<<"Error: pvalueArray and labelArray have different number of tuples.\n";
          abort();
        }


      double n1 = 0;
      double n2 = 0;

      for(int i = 0; i < n; i++)
        {
          if (fabs(labelArray->GetValue(i) - v1) < 1e-3)
            {
              ++n1;
              if ( pvalueArray->GetValue(i) <= 0.04)
                {
                  ++n2;
                }
            }
        }

      std::cout << "ratio is: " << n2/n1 << std::endl;
    }



  return EXIT_SUCCESS;
}
