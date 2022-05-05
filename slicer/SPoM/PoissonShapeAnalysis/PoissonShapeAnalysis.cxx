#include "lib/utilitiesIO.h"
#include "lib/SumOfTwoPoissonShapeAnalysis.h"

#include "itkPluginUtilities.h"

#include "PoissonShapeAnalysisCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

// template <typename TPixel>
// int DoIt( int argc, char * argv[], TPixel )
// {
//   PARSE_ARGS;

//   typedef TPixel InputPixelType;
//   typedef TPixel OutputPixelType;

//   const unsigned int Dimension = 3;

//   typedef itk::Image<InputPixelType,  Dimension> InputImageType;
//   typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

//   typedef itk::ImageFileReader<InputImageType>  ReaderType;

//   typename ReaderType::Pointer reader = ReaderType::New();

//   reader->SetFileName( inputVolume.c_str() );

//   typedef itk::SmoothingRecursiveGaussianImageFilter<
//     InputImageType, OutputImageType>  FilterType;
//   typename FilterType::Pointer filter = FilterType::New();
//   filter->SetInput( reader->GetOutput() );
//   filter->SetSigma( sigma );

//   typedef itk::ImageFileWriter<OutputImageType> WriterType;
//   typename WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( outputVolume.c_str() );
//   writer->SetInput( filter->GetOutput() );
//   writer->SetUseCompression(1);
//   writer->Update();

//   return EXIT_SUCCESS;
// }

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // if (argc < 5)
  //   {
  //     std::cerr<<"args: shapeList1 shapeList2 meanImageName meanVTPName [doRegistration 0/1] [meanShape]\n";
  //     exit(-1);
  //   }

  // std::string shapeList1Name(argv[1]);
  // std::string shapeList2Name(argv[2]);
  // std::string meanImageName(argv[3]);
  // std::string meanVTKName(argv[4]);

  // bool doRegistration = true;
  // if (argc > 5)
  //   {
  //     if (atoi(argv[5]) == 0)
  //       {
  //         doRegistration = false;
  //       }
  //   }

  const unsigned int ImageDim = 3;
  typedef short pixel_t;
  typedef itk::Image<pixel_t, ImageDim> image_t;

  // Read in training images' names
  // std::vector< std::string > shapeNameList1 = ShapeAnalysis::readTextLineToListOfString<char>(shapeList1Name.c_str());
  // std::vector< std::string > shapeNameList2 = ShapeAnalysis::readTextLineToListOfString<char>(shapeList2Name.c_str());

  std::vector< image_t::Pointer > shapeList1 = ShapeAnalysis::readImageSeries<image_t>(shapeNameList1);
  std::vector< image_t::Pointer > shapeList2 = ShapeAnalysis::readImageSeries<image_t>(shapeNameList2);

  typedef ShapeAnalysis::SumOfTwoPoissonShapeAnalysis<image_t> SumOfTwoPoissonShapeAnalysisType;
  SumOfTwoPoissonShapeAnalysisType foo;

  foo.SetShapeGroup1(shapeList1);
  foo.SetShapeGroup2(shapeList2);
  foo.setPerformRegistration(doRegistration);

  // if (argc > 6)
  //   {
  //     foo.setMeanShapeImage(ShapeAnalysis::readImage<SumOfTwoPoissonShapeAnalysisType::FloatImageType>(argv[6]));
  //   }

  foo.update();

  ShapeAnalysis::writeImage<SumOfTwoPoissonShapeAnalysisType::FloatImageType>(foo.getMeanShapeImage(), meanImageName.c_str());
  ShapeAnalysis::writePolydataToXMLFile<char>(meanVTKName.c_str(), foo.getMeanShapeSurface());

  return EXIT_SUCCESS;
}


// {
//   PARSE_ARGS;

//   itk::ImageIOBase::IOPixelType     pixelType;
//   itk::ImageIOBase::IOComponentType componentType;

//   try
//     {
//     itk::GetImageType(inputVolume, pixelType, componentType);

//     // This filter handles all types on input, but only produces
//     // signed types
//     switch( componentType )
//       {
//       case itk::ImageIOBase::UCHAR:
//         return DoIt( argc, argv, static_cast<unsigned char>(0) );
//         break;
//       case itk::ImageIOBase::CHAR:
//         return DoIt( argc, argv, static_cast<signed char>(0) );
//         break;
//       case itk::ImageIOBase::USHORT:
//         return DoIt( argc, argv, static_cast<unsigned short>(0) );
//         break;
//       case itk::ImageIOBase::SHORT:
//         return DoIt( argc, argv, static_cast<short>(0) );
//         break;
//       case itk::ImageIOBase::UINT:
//         return DoIt( argc, argv, static_cast<unsigned int>(0) );
//         break;
//       case itk::ImageIOBase::INT:
//         return DoIt( argc, argv, static_cast<int>(0) );
//         break;
//       case itk::ImageIOBase::ULONG:
//         return DoIt( argc, argv, static_cast<unsigned long>(0) );
//         break;
//       case itk::ImageIOBase::LONG:
//         return DoIt( argc, argv, static_cast<long>(0) );
//         break;
//       case itk::ImageIOBase::FLOAT:
//         return DoIt( argc, argv, static_cast<float>(0) );
//         break;
//       case itk::ImageIOBase::DOUBLE:
//         return DoIt( argc, argv, static_cast<double>(0) );
//         break;
//       case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
//       default:
//         std::cerr << "Unknown input image pixel component type: ";
//         std::cerr << itk::ImageIOBase::GetComponentTypeAsString( componentType );
//         std::cerr << std::endl;
//         return EXIT_FAILURE;
//         break;
//       }
//     }

//   catch( itk::ExceptionObject & excep )
//     {
//     std::cerr << argv[0] << ": exception caught !" << std::endl;
//     std::cerr << excep << std::endl;
//     return EXIT_FAILURE;
//     }
//   return EXIT_SUCCESS;
// }
