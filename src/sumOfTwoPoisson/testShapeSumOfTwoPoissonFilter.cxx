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

#include "ShapeSumOfTwoPoissonFilter.h"

int main(int argc, char** argv)
{
  if (argc < 3)
    {
      std::cerr<<"args: inputShape sumOfTwoPoisson\n";
      exit(-1);
    }

  std::string inputShapeName(argv[1]);
  std::string sumOfTwoPoissonName(argv[2]);

  const unsigned int ImageDim = 3;
  typedef short pixel_t;
  typedef itk::Image<pixel_t, ImageDim> image_t;

  // Read in training images' names
  image_t::Pointer inputShape = ShapeAnalysis::readImage<image_t>(inputShapeName.c_str());

  typedef ShapeAnalysis::ShapeSumOfTwoPoissonFilter<image_t> ShapeSumOfTwoPoissonFilterType;
  ShapeSumOfTwoPoissonFilterType foo;

  foo.setInputShape(inputShape);
  foo.update();

  ShapeSumOfTwoPoissonFilterType::FloatImagePointer sumOfPoisson = foo.getPoissonImage();

  ShapeAnalysis::writeImage<ShapeSumOfTwoPoissonFilterType::FloatImageType>(sumOfPoisson, sumOfTwoPoissonName.c_str());


  return 0;
}
