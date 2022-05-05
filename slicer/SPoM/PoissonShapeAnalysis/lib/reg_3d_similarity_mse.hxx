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


#ifndef reg_3d_similarity_mse_hxx_
#define reg_3d_similarity_mse_hxx_

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


// itk
#include "itkSimilarity3DTransform.h"

#include "itkCenteredTransformInitializer.h"

#include "itkCommand.h"

#include "itkImage.h"
#include "itkImageRegistrationMethod.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"


// local
#include "reg_3d_similarity_mse.h"
// #include "utilitiesIO.h"
// #include "utilitiesImage.h"

namespace ShapeAnalysis
{
  /**********************************************************************************
   *
   *********************************************************************************/
  template<typename fix_image_t, typename moving_image_t>
  itk::Similarity3DTransform<double>::Pointer
  similarityMSERegistration(typename fix_image_t::Pointer fixImg, typename moving_image_t::Pointer movingImg, double& finalCost, int numThreads)
  {
    typedef itk::Similarity3DTransform< double > transform_t;

    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::MeanSquaresImageToImageMetric< fix_image_t, moving_image_t >  MetricType;
    typedef itk::LinearInterpolateImageFunction< moving_image_t, double >    InterpolatorType;
    typedef itk::ImageRegistrationMethod< fix_image_t, moving_image_t >    RegistrationType;

    typename MetricType::Pointer         metric        = MetricType::New();
    typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
    typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    typename RegistrationType::Pointer   registration  = RegistrationType::New();

    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );

    if (numThreads != -1 && numThreads >= 1 && numThreads <= 32)
      {
        // I don't have SMP with 32 cores.
        registration->SetNumberOfThreads(numThreads);
      }

    typename transform_t::Pointer  transform = transform_t::New();
    registration->SetTransform( transform );

    registration->SetFixedImage(  fixImg    );
    registration->SetMovingImage( movingImg );

    registration->SetFixedImageRegion( fixImg->GetBufferedRegion() );

    typedef itk::CenteredTransformInitializer< transform_t, fix_image_t, moving_image_t >  TransformInitializerType;
    typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransform(   transform );
    initializer->SetFixedImage(  fixImg );
    initializer->SetMovingImage( movingImg );
    initializer->MomentsOn();
    //initializer->GeometryOn();
    initializer->InitializeTransform();


    registration->SetInitialTransformParameters( transform->GetParameters() );

    double translationScale = 1.0 / 1000.0;

    typedef OptimizerType::ScalesType       OptimizerScalesType;
    OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

    optimizerScales[0] = 1.0; // versor.x
    optimizerScales[1] = 1.0; // versor.y
    optimizerScales[2] = 1.0; // versor.z
    optimizerScales[3] = translationScale; // translation.x
    optimizerScales[4] = translationScale; // translation.y
    optimizerScales[5] = translationScale; // translation.z
    optimizerScales[6] = 1.0; // scale

    optimizer->SetScales( optimizerScales );

    double steplength = 0.1;
    optimizer->SetMaximumStepLength( steplength ); 
    optimizer->SetMinimumStepLength( 0.0005 );

    unsigned int maxNumberOfIterations = 2000;
    optimizer->SetNumberOfIterations( maxNumberOfIterations );
    optimizer->MinimizeOn();

    try 
      { 
        //registration->StartRegistration(); 
        registration->Update(); 
        //     std::cout << "Optimizer stop condition: "
        //               << registration->GetOptimizer()->GetStopConditionDescription()
        //               << std::endl;
      } 
    catch( itk::ExceptionObject & err ) 
      { 
        std::cerr << "ExceptionObject caught !" << std::endl; 
        std::cerr << err << std::endl; 
        exit(-1);
      } 

    // record final cost function value
    finalCost = optimizer->GetValue();
  
    //   //tst
    //   std::cout<<"finalCostValue = "<<finalCostValue<<std::endl;
    //   //tst//

    OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

    typename transform_t::Pointer finalTransform = transform_t::New();
    finalTransform->SetParameters( finalParameters );
    finalTransform->SetFixedParameters( transform->GetFixedParameters() );

    return finalTransform;
  }

} // namespace ShapeAnalysis

#endif
