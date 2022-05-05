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


#ifndef SumOfTwoPoissonShapeAnalysis_hxx_
#define SumOfTwoPoissonShapeAnalysis_hxx_

#include <vector>

// alglib
#include "ap.h"
#include "statistics.h"

// itk
#include "itkAddImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkPoint.h"
#include "itkResampleImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"


// itkVTKGlue
#include "itkImageToVTKImageFilter.h"

// vtk
#include "vtkDoubleArray.h"
#include "vtkMarchingCubes.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkTriangleFilter.h"
#include "vtkButterflySubdivisionFilter.h"
#include "vtkSmoothPolyDataFilter.h"

// local
#include "reg_3d_similarity_mse.h"
#include "transformImage.h"
#include "utilitiesImage.h"
#include "ShapeSumOfTwoPoissonFilter.h"

#include "SumOfTwoPoissonShapeAnalysis.h"

#include "itkVtkMeshConversion.h"


#include "utilitiesImage.h"

// for debug
#include "utilitiesIO.h"
#include "utilitiesVTK.h"

namespace ShapeAnalysis
{
  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::update()
  {
    _preprocess();

    _computeMeanShape();

    _computeMeanShapeSurface();

    _computePoissonMaps();

    _poissonDistanceMapOnMeanShapeSurface();

    _computePValueOnMeanShape();

    m_allDone = true;

    return;
  }

  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_preprocess()
  {
    _cropInputImages();

    // //--------------------------------------------------------------------------------
    // // Check image regions
    // if (m_useOutsideMeanShape)
    //   {
    //     if (m_meanShape->GetLargestPossibleRegion() != m_shapeGroup1[0]->GetLargestPossibleRegion())
    //       {
    //         std::cerr<<"Error: input mean shape does not match region of the input shape.\n";
    //         std::cerr<<"m_meanShape->GetLargestPossibleRegion():\n";
    //         std::cerr<<m_meanShape->GetLargestPossibleRegion();

    //         std::cerr<<"================================================================================\n";
    //         std::cerr<<"m_shapeGroup1[0]->GetLargestPossibleRegion():\n";
    //         std::cerr<<m_shapeGroup1[0]->GetLargestPossibleRegion();

    //         abort();
    //       }
    //   }
    // //================================================================================



    /**
     * The input shapes, binary masks of general segmentations, are
     * not isotropic. To avoid possible problem of itk-vtk image
     * transfer, marching cube, and sample on the mean shape surface,
     * let me first isotropicize the images. To that end, we just need
     * to do that for the first image because the registration will
     * make the others the same.
     */
    //_isotropicizeFirstImage();

    return;
  }

  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_cropInputImages()
  {
    //dbg
    char group1Name[1000];
    char group2Name[1000];
    //dbg, end


    for (std::size_t it = 0; it < m_shapeGroup1.size(); ++it)
      {
        m_shapeGroup1[it] = cropNonZeroRegionFromImage<ShapeImageType>(m_shapeGroup1[it]);

    //dbg
        sprintf(group1Name, "group1_%lu.nrrd", it);
        writeImage<ShapeImageType>(m_shapeGroup1[it], group1Name);
    //dbg, end
      }

    for (std::size_t it = 0; it < m_shapeGroup2.size(); ++it)
      {
        m_shapeGroup2[it] = cropNonZeroRegionFromImage<ShapeImageType>(m_shapeGroup2[it]);

    //dbg
        sprintf(group2Name, "group2_%lu.nrrd", it);
        writeImage<ShapeImageType>(m_shapeGroup2[it], group2Name);
    //dbg, end
      }

    return;
  }



  template< typename TShapeImageType >
  void
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_isotropicizeFirstImage()
  {
    if (m_shapeGroup1.empty())
      {
        std::cerr<<"Error: m_shapeGroup1.empty()\n";
        abort();
      }


    //dbg
    writeImage<ShapeImageType>(m_shapeGroup1[0], "firstImageBefore.nrrd");
    //dbg, end

    typedef typename ShapeImageType::SpacingType SpacingType;
    SpacingType inputSpacing = m_shapeGroup1[0]->GetSpacing();
    float smallestSpacing = inputSpacing[0] > inputSpacing[1]?inputSpacing[1]:inputSpacing[0];
    smallestSpacing = smallestSpacing > inputSpacing[2]?inputSpacing[2]:smallestSpacing;

    SpacingType outputSpacing;
    outputSpacing.Fill(smallestSpacing);

    typedef typename ShapeImageType::SizeType SizeType;
    SizeType inputSize = m_shapeGroup1[0]->GetLargestPossibleRegion().GetSize();
    SizeType outputSize;
    outputSize[0] = static_cast<typename SizeType::SizeValueType>(inputSize[0]*(inputSpacing[0]/outputSpacing[0]));
    outputSize[1] = static_cast<typename SizeType::SizeValueType>(inputSize[1]*(inputSpacing[1]/outputSpacing[1]));
    outputSize[2] = static_cast<typename SizeType::SizeValueType>(inputSize[2]*(inputSpacing[2]/outputSpacing[2]));

    //dbg
    std::cout<<inputSize<<"\t"<<outputSize<<std::endl;
    //dbg, end

    typedef itk::IdentityTransform<double, Dim> IDTransformType; // has to be double, float causes compile error, not sure why
    typedef itk::ResampleImageFilter<ShapeImageType, ShapeImageType> ResampleImageFilterType;
    typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
    resample->SetInput(m_shapeGroup1[0]);
    resample->SetOutputOrigin(m_shapeGroup1[0]->GetOrigin());
    resample->SetSize(outputSize);
    resample->SetOutputSpacing(outputSpacing);
    resample->SetTransform(IDTransformType::New());
    resample->UpdateLargestPossibleRegion();

    m_shapeGroup1[0] = resample->GetOutput();

    //dbg
    writeImage<ShapeImageType>(m_shapeGroup1[0], "firstImageAfter.nrrd");
    exit(0);
    //dbg, end

    return;
  }


  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_computeMeanShape()
  {
    /**
     * Given a lsit of binary shapes, register them to the first one,
     * using similarity transformation and MSE metric, to get the mean
     * shape. Since the similarity trans does NOT change the shape,
     * then, register to any shape does not change the shape of the mean
     * shape. Then, register eash shape, again, to the mean shape, and
     * compute mean shape again. While, that may be enough, let's not
     * register to the new mean shape again and again till convergence.
     */
    if (m_useOutsideMeanShape)
      {
        m_registrationFixedImage = m_meanShape;
      }
    else
      {
        m_registrationFixedImage = castItkImage<typename ShapeImageType::PixelType, typename FloatImageType::PixelType>( m_shapeGroup1[0] );

        m_meanShape = FloatImageType::New();
        m_meanShape->SetRegions(m_shapeGroup1[0]->GetLargestPossibleRegion() );
        m_meanShape->Allocate();
        m_meanShape->FillBuffer(0.0);
        m_meanShape->CopyInformation(m_shapeGroup1[0]);
      }

    // FloatImagePointer tmpMeanShape = FloatImageType::New();
    // tmpMeanShape->SetRegions(m_shapeGroup1[0]->GetLargestPossibleRegion() );
    // tmpMeanShape->Allocate();
    // tmpMeanShape->FillBuffer(0.0);
    // tmpMeanShape->CopyInformation(m_shapeGroup1[0]);

    long n1 = static_cast<long>(m_shapeGroup1.size());
    long n2 = static_cast<long>(m_shapeGroup2.size());

    long n = n1 + n2;

    std::cout<<"similarity registration of group-1 (totally "<<n1<<" ) shapes.\n   "<<std::flush;
    for (long it = 0; it < n1; ++it)
      {
        std::cout<<it<<", "<<std::flush;

        double finalCost = 0.0;
        typedef itk::Similarity3DTransform<double> SimilarityTransformType;

        typename SimilarityTransformType::Pointer trans = SimilarityTransformType::New();
        if (m_performRegistration)
          {
            trans = similarityMSERegistration<FloatImageType, ShapeImageType>(m_registrationFixedImage, m_shapeGroup1[it], finalCost);
          }
        else
          {
            trans->SetIdentity();
          }

        FloatType fillInValue = 0.0;
        char interpolationType = 1; // linear interp

        FloatImagePointer registeredShape                               \
          = transformImage<ShapeImageType, FloatImageType, FloatImageType>(trans, m_shapeGroup1[it], m_registrationFixedImage, fillInValue, interpolationType);

        m_similarityRegisteredShapeGroup1.push_back(registeredShape);

        if (!m_useOutsideMeanShape)
          {
            typedef itk::AddImageFilter< FloatImageType, FloatImageType, FloatImageType > AddImageFilterType;
            typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
            addFilter->SetInput1(m_meanShape);
            addFilter->SetInput2(registeredShape);
            addFilter->Update();

            m_meanShape = addFilter->GetOutput();
          }
      }
    std::cout<<std::endl;

    std::cout<<"similarity registration of group-2 (totally "<<n2<<" ) shapes.\n   "<<std::flush;
    for (long it = 0; it < n2; ++it)
      {
        std::cout<<it<<", "<<std::flush;

        double finalCost = 0.0;
        typedef itk::Similarity3DTransform<double> SimilarityTransformType;

        typename SimilarityTransformType::Pointer trans = SimilarityTransformType::New();
        if (m_performRegistration)
          {
            trans = similarityMSERegistration<FloatImageType, ShapeImageType>(m_registrationFixedImage, m_shapeGroup2[it], finalCost);
          }
        else
          {
            trans->SetIdentity();
          }

        FloatType fillInValue = 0.0;
        char interpolationType = 1; // linear interp

        FloatImagePointer registeredShape                               \
          = transformImage<ShapeImageType, FloatImageType, FloatImageType>(trans, m_shapeGroup2[it], m_registrationFixedImage, fillInValue, interpolationType);

        m_similarityRegisteredShapeGroup2.push_back(registeredShape);

        if (!m_useOutsideMeanShape)
          {
            typedef itk::AddImageFilter< FloatImageType, FloatImageType, FloatImageType > AddImageFilterType;
            typename AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
            addFilter->SetInput1(m_meanShape);
            addFilter->SetInput2(registeredShape);
            addFilter->Update();

            m_meanShape = addFilter->GetOutput();
          }
      }
    std::cout<<std::endl;

    // divide meanShape by n, the divideByConstantImageFilter seems to be in the deprecated category? so i just do it
    if (!m_useOutsideMeanShape)
      {
        typedef itk::ImageRegionIterator<FloatImageType> ImageRegionIterator;

        FloatType nn = static_cast<FloatType>(n);

        ImageRegionIterator it(m_meanShape, m_meanShape->GetLargestPossibleRegion());
        for (it.GoToBegin(); !it.IsAtEnd(); ++it)
          {
            it.Set(it.Get()/nn);
          }
      }


    /**
     * Now, we have the mean shape and all the shapes are registered
     * to it with similarity trans.
     */

    return;
  }


  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_computeMeanShapeSurface()
  {
    typename FloatImageType::PixelType lowerThreshold = -1e5;
    typename FloatImageType::PixelType upperThreshold = 0.5;

    typedef itk::BinaryThresholdImageFilter< FloatImageType, FloatImageType > BinaryThresholdFilterType;
    typename BinaryThresholdFilterType::Pointer threshold = BinaryThresholdFilterType::New();
    threshold->SetInput( m_meanShape );
    threshold->SetLowerThreshold( lowerThreshold );
    threshold->SetUpperThreshold( upperThreshold );
    threshold->SetOutsideValue( 1 );
    threshold->Update();

    //typedef itk::Mesh< double, FloatImageType::ImageDimension > MeshType;

    typedef itk::BinaryMask3DMeshSource< FloatImageType, MeshType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( threshold->GetOutput() );
    filter->SetObjectValue( 1 );
    filter->Update();

    vtkSmartPointer<vtkPolyData> pd = ITKMeshToVtkPolyData( filter->GetOutput() );

    // // //dbg
    // writePolydataToXMLFile<char>("pd.vtp", pd);
    // // //dbg, end

    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangleFilter->SetInputData(pd);
    triangleFilter->Update();

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(triangleFilter->GetOutput());
    smoothFilter->SetNumberOfIterations(300);
    smoothFilter->Update();

    int numberOfSubdivisions = 2;
    vtkSmartPointer<vtkButterflySubdivisionFilter> butterflySubdivisionFilter = vtkSmartPointer<vtkButterflySubdivisionFilter>::New();
    butterflySubdivisionFilter->SetNumberOfSubdivisions(numberOfSubdivisions);
    //butterflySubdivisionFilter->SetInputData(reverseSense->GetOutput());
    butterflySubdivisionFilter->SetInputData(smoothFilter->GetOutput());
    butterflySubdivisionFilter->Update();


    m_meanShapeSurface = butterflySubdivisionFilter->GetOutput();

    // //dbg
    writePolydataToXMLFile<char>("m_meanShapeSurface.vtp", m_meanShapeSurface);
    // std::cout<<mc->GetOutput()->GetNumberOfPoints()<<std::endl;
    // std::cout<<mc->GetOutput()->GetNumberOfCells()<<std::endl;
    // //dbg, end

    return;
  }



  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_computePoissonMaps()
  {
    typedef itk::Image<char, Dim> CharImageType;

    std::cout<<"Computeing Poisson map of group 1, totally "<<m_similarityRegisteredShapeGroup1.size()<<" shapes. \n   "<<std::flush;
    m_poissonMapGroup1.resize(m_similarityRegisteredShapeGroup1.size());

    for (std::size_t it = 0; it < m_similarityRegisteredShapeGroup1.size(); ++it)
      {
        std::cout<<it<<", "<<std::flush;

        // //dbg
        // char fileNameB[255];
        // sprintf(fileNameB, "bin1_%ld.nrrd", it);
        // writeImage<CharImageType>(binarilizeImage<FloatImageType, CharImageType>(m_similarityRegisteredShapeGroup1[it], 0.5, 1e10, 1, 0), fileNameB);
        // //dbg, end

        typedef ShapeAnalysis::ShapeSumOfTwoPoissonFilter<CharImageType> ShapeSumOfTwoPoissonFilterType;
        ShapeSumOfTwoPoissonFilterType foo;

        foo.setInputShape(binarilizeImage<FloatImageType, CharImageType>(m_similarityRegisteredShapeGroup1[it], 0.5, 1e10, 1, 0));
        foo.update();

        m_poissonMapGroup1[it] = foo.getPoissonImage();

        // //dbg
        // char fileName[255];
        // sprintf(fileName, "poisson1_%ld.nrrd", it);
        // writeImage<FloatImageType>(m_poissonMapGroup1[it], fileName);
        // //dbg, end
      }
    std::cout<<std::endl;

    std::cout<<"Computeing Poisson map of group 2, totally "<<m_similarityRegisteredShapeGroup2.size()<<" shapes. \n   "<<std::flush;
    m_poissonMapGroup2.resize(m_similarityRegisteredShapeGroup2.size());

    for (std::size_t it = 0; it < m_similarityRegisteredShapeGroup2.size(); ++it)
      {
        std::cout<<it<<", "<<std::flush;

        typedef ShapeAnalysis::ShapeSumOfTwoPoissonFilter<CharImageType> ShapeSumOfTwoPoissonFilterType;
        ShapeSumOfTwoPoissonFilterType foo;

        foo.setInputShape(binarilizeImage<FloatImageType, CharImageType>(m_similarityRegisteredShapeGroup2[it], 0.5, 1e10, 1, 0));
        foo.update();

        // //dbg
        // char fileNameB[255];
        // sprintf(fileNameB, "bin2_%ld.nrrd", it);
        // writeImage<CharImageType>(binarilizeImage<FloatImageType, CharImageType>(m_similarityRegisteredShapeGroup2[it], 0.5, 1e10, 1, 0), fileNameB);
        // //dbg, end

        m_poissonMapGroup2[it] = foo.getPoissonImage();

        // //dbg
        // char fileName[255];
        // sprintf(fileName, "poisson2_%ld.nrrd", it);
        // writeImage<FloatImageType>(m_poissonMapGroup2[it], fileName);
        // //dbg, end
      }

    std::cout<<std::endl;

    return;
  }

  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_poissonDistanceMapOnMeanShapeSurface()
  {
    /**
     * Put the mean shape surface to each det jacobian field, and get
     * values on mean shape surface.
     */

    // Write all of the coordinates of the points in the vtkPolyData to the console.
    vtkIdType n = m_meanShapeSurface->GetNumberOfPoints();

    // //dbg
    // std::cout<<"number of points on mean surface = "<<n<<std::endl;
    // //dbg, end

    typedef itk::Point<FloatType, Dim> ITKPointType;
    typedef itk::ContinuousIndex<FloatType, Dim> ContinuousIndexType;

    typedef itk::LinearInterpolateImageFunction<FloatImageType, FloatType> LinearInterpolateImageFunctionType;

    for (std::size_t it = 0; it < m_poissonMapGroup1.size(); ++it)
      {
        typename LinearInterpolateImageFunctionType::Pointer interpolator = LinearInterpolateImageFunctionType::New();
        interpolator->SetInputImage(m_poissonMapGroup1[it]);

        std::vector<FloatType> v(n, 0.0);

        for(vtkIdType i = 0; i < n; ++i)
          {
            double p[3];
            m_meanShapeSurface->GetPoint(i, p);

            ITKPointType pos;
            pos[0] = p[0];
            pos[1] = p[1];
            pos[2] = p[2];

            ContinuousIndexType cidx;
            m_poissonMapGroup1[it]->TransformPhysicalPointToContinuousIndex(pos, cidx);

            if (interpolator->IsInsideBuffer(cidx))
              {
                v[i] = interpolator->EvaluateAtContinuousIndex(cidx);
              }
            else
              {
                v[i] = -100.0;
              }

            // dbg
            //std::cout<<pos<<"\t"<<cidx<<std::endl;
            // dbg, end

            // // dbg
            // v[i] = 1.0 + interpolator->EvaluateAtContinuousIndex(cidx);
            // // dbg, end

            // // dbg
            // std::cout<<v[i]<<" ";
            // // dbg, pause
          }

        // // dbg, cont
        // std::cout<<std::endl;
        // // dbg, end

        // //dbg
        // appendScalarFieldToPolydataPoints<FloatType>(m_meanShapeSurface, v);
        // char fileName[255];
        // sprintf(fileName, "group1_%ld.vtp", it);
        // writePolydataToXMLFile<char>(fileName, m_meanShapeSurface);
        // //dbg, end

        m_poissonDistanceOnMeanSurface1.push_back(v);
      }


    for (std::size_t it = 0; it < m_poissonMapGroup2.size(); ++it)
      {
        typename LinearInterpolateImageFunctionType::Pointer interpolator = LinearInterpolateImageFunctionType::New();
        interpolator->SetInputImage(m_poissonMapGroup2[it]);

        std::vector<FloatType> v(n, 0.0);

        for(vtkIdType i = 0; i < n; ++i)
          {
            double p[3];
            m_meanShapeSurface->GetPoint(i, p);

            ITKPointType pos;
            pos[0] = p[0];
            pos[1] = p[1];
            pos[2] = p[2];

            ContinuousIndexType cidx;
            m_poissonMapGroup2[it]->TransformPhysicalPointToContinuousIndex(pos, cidx);

            if (interpolator->IsInsideBuffer(cidx))
              {
                v[i] = interpolator->EvaluateAtContinuousIndex(cidx);
              }
            else
              {
                v[i] = -100.0;
              }

            // dbg
            //std::cout<<pos<<"\t"<<cidx<<std::endl;
            // dbg, end

            // // dbg
            // v[i] = 1.0 + interpolator->EvaluateAtContinuousIndex(cidx);
            // // dbg, end

            // // dbg
            // std::cout<<v[i]<<" ";
            // // dbg, pause
          }

        // // dbg, cont
        // std::cout<<std::endl;
        // // dbg, end

        // //dbg
        // appendScalarFieldToPolydataPoints<FloatType>(m_meanShapeSurface, v);
        // char fileName[255];
        // sprintf(fileName, "group2_%ld.vtp", it);
        // writePolydataToXMLFile<char>(fileName, m_meanShapeSurface);
        // //dbg, end

        m_poissonDistanceOnMeanSurface2.push_back(v);
      }




    // //dbg
    // ContinuousIndexType tmpCIdx;
    // tmpCIdx[0] = 0;
    // tmpCIdx[1] = 0;
    // tmpCIdx[2] = 0;

    // ITKPointType tmpPos;

    // //m_poissonMapGroup2[0]->TransformPhysicalPointToContinuousIndex(pos, cidx);
    // m_poissonMapGroup2[0]->TransformPhysicalPointToContinuousIndex(pos, cidx);

    //dbg, end


    return;
  }


  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_poissonDistance(FloatImagePointer poissonField, double startingPosition[3])
  {
    // poissonField should NOT be edited here, infact it should be
    // const pointer. startingPosition is in physical space.

    // TODO, add the impl here. for Poisson, the difficulty in the
    // "monotonicity" is not solved. discuss

  }

  template< typename TShapeImageType >
  void SumOfTwoPoissonShapeAnalysis<TShapeImageType>::_computePValueOnMeanShape()
  {
    vtkIdType n = m_meanShapeSurface->GetNumberOfPoints();
    m_pValuesOnMeanShapeSurface.resize(n);


    alglib::ae_int_t n1 = static_cast<alglib::ae_int_t>(m_poissonDistanceOnMeanSurface1.size());
    alglib::ae_int_t n2 = static_cast<alglib::ae_int_t>(m_poissonDistanceOnMeanSurface2.size());

    alglib::real_1d_array x;
    alglib::real_1d_array y;

    double* xx = new double[n1];
    double* yy = new double[n2];

    double bothtails = 0;
    double lefttail = 0;
    double righttail = 0;

    for(vtkIdType i = 0; i < n; ++i)
      {
        for (alglib::ae_int_t i1 = 0; i1 < n1; ++i1)
          {
            xx[i1] = m_poissonDistanceOnMeanSurface1[i1][i];
          }

        for (alglib::ae_int_t i2 = 0; i2 < n2; ++i2)
          {
            yy[i2] = m_poissonDistanceOnMeanSurface2[i2][i];
          }

        x.setcontent(n1, xx);
        y.setcontent(n2, yy);

        alglib::studentttest2(x, n1, y, n2, bothtails, lefttail, righttail);

        m_pValuesOnMeanShapeSurface[i] = bothtails;
      }

    delete[] xx;
    delete[] yy;

    // put p-values on the mean shape surface
    vtkSmartPointer<vtkDoubleArray> pvalues = vtkSmartPointer<vtkDoubleArray>::New();
    pvalues->SetNumberOfValues(n);

    for(vtkIdType i = 0; i < n; ++i)
      {
        pvalues->SetValue(i, m_pValuesOnMeanShapeSurface[i]);
      }
    pvalues->SetName("pValue");
    m_meanShapeSurface->GetPointData()->SetScalars(pvalues);

    return;
  }


  template< typename TShapeImageType >
  void
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::setMeanShapeImage(FloatImagePointer meanshape)
  {
    m_meanShape = meanshape;

    m_useOutsideMeanShape = true;
  }


  template< typename TShapeImageType >
  typename SumOfTwoPoissonShapeAnalysis<TShapeImageType>::FloatImagePointer
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::getMeanShapeImage()
  {
    if (m_useOutsideMeanShape || m_allDone)
      {
        return m_meanShape;
      }
    else
      {
        update();

        return getMeanShapeImage();
      }
  }

  template< typename TShapeImageType >
  vtkSmartPointer<vtkPolyData>
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::getMeanShapeSurface()
  {
    /* Even if mean shape is from outside, the mean shape surface has
       to be computed. */
    if (m_allDone)
      {
        return m_meanShapeSurface;
      }
    else
      {
        update();
        return getMeanShapeSurface();
      }
  }


  template< typename TShapeImageType >
  void
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::SetShapeGroup1(const std::vector< ShapeImageTypePointer >& shapeGroup)
  {
    m_shapeGroup1 = shapeGroup;

    // TODO check size is 0?

    return;
  }

  template< typename TShapeImageType >
  void
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::SetShapeGroup2(const std::vector< ShapeImageTypePointer >& shapeGroup)
  {
    m_shapeGroup2 = shapeGroup;

    // TODO check size is 0?

    return;
  }


  template< typename TShapeImageType >
  void
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::setPerformRegistration(bool reg)
  {
    m_performRegistration = reg;

    return;
  }


  template< typename TShapeImageType >
  SumOfTwoPoissonShapeAnalysis<TShapeImageType>::SumOfTwoPoissonShapeAnalysis()
  {
    m_allDone = false;

    m_performRegistration = true;

    m_meanShape = 0;
    m_useOutsideMeanShape = false;

    return;
  }





}// ShapeAnalysis

#endif
