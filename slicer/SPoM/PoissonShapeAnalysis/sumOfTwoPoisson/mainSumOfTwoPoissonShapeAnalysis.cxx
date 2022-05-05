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

#include "SumOfTwoPoissonShapeAnalysis.h"

int main(int argc, char** argv)
{
  if (argc < 5)
    {
      std::cerr<<"args: shapeList1 shapeList2 meanImageName meanVTPName [doRegistration 0/1] [meanShape]\n";
      exit(-1);
    }

  std::string shapeList1Name(argv[1]);
  std::string shapeList2Name(argv[2]);
  std::string meanImageName(argv[3]);
  std::string meanVTKName(argv[4]);

  bool doRegistration = true;
  if (argc > 5)
    {
      if (atoi(argv[5]) == 0)
        {
          doRegistration = false;
        }
    }

  const unsigned int ImageDim = 3;
  typedef short pixel_t;
  typedef itk::Image<pixel_t, ImageDim> image_t;

  // Read in training images' names
  std::vector< std::string > shapeNameList1 = ShapeAnalysis::readTextLineToListOfString<char>(shapeList1Name.c_str());
  std::vector< std::string > shapeNameList2 = ShapeAnalysis::readTextLineToListOfString<char>(shapeList2Name.c_str());

  std::vector< image_t::Pointer > shapeList1 = ShapeAnalysis::readImageSeries<image_t>(shapeNameList1);
  std::vector< image_t::Pointer > shapeList2 = ShapeAnalysis::readImageSeries<image_t>(shapeNameList2);

  typedef ShapeAnalysis::SumOfTwoPoissonShapeAnalysis<image_t> SumOfTwoPoissonShapeAnalysisType;
  SumOfTwoPoissonShapeAnalysisType foo;

  foo.SetShapeGroup1(shapeList1);
  foo.SetShapeGroup2(shapeList2);
  foo.setPerformRegistration(doRegistration);

  if (argc > 6)
    {
      foo.setMeanShapeImage(ShapeAnalysis::readImage<SumOfTwoPoissonShapeAnalysisType::FloatImageType>(argv[6]));
    }

  foo.update();

  ShapeAnalysis::writeImage<SumOfTwoPoissonShapeAnalysisType::FloatImageType>(foo.getMeanShapeImage(), meanImageName.c_str());
  ShapeAnalysis::writePolydataToXMLFile<char>(meanVTKName.c_str(), foo.getMeanShapeSurface());

  return 0;
}
