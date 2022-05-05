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


#ifndef reg_3d_rigid_mse_h_
#define reg_3d_rigid_mse_h_

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkVersorRigid3DTransform.h"

namespace ShapeAnalysis
{
  template<typename fix_image_t, typename moving_image_t>
  itk::VersorRigid3DTransform<double>::Pointer
  rigidMSERegistration(typename fix_image_t::Pointer fixedImage, typename moving_image_t::Pointer movingImage, double& finalCost, int numThreads = -1);

} //ShapeAnalysis

#include "reg_3d_rigid_mse.hxx"

#endif
