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


#ifndef utilitiesIO_hxx_
#define utilitiesIO_hxx_

#include <csignal>
#include <string>

// itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"

#include "itkVector.h"

#include "itkTransformFactoryBase.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"



// vnl
#include "vnl/vnl_matrix.h"

// vtk
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"

#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

// local
#include "utilitiesIO.h"

namespace ShapeAnalysis
{
  /**
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName)
  {
    typedef itk::ImageFileReader< itkImage_t > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename itkImage_t::Pointer image;

    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl;
        std::cerr<< err << std::endl;
        raise(SIGABRT);
      }

    return image;
  }


  /**
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName)
  {
    typedef itk::ImageFileWriter< itkImage_t > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }
  }


  /**
   * Read a series of images.
   */
  template< typename itkImage_t >
  std::vector< typename itkImage_t::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList )
  {
    typedef typename itkImage_t::Pointer itkImagePointer_t;
    typedef std::vector< itkImagePointer_t > itkImageList_t;
    typedef itk::ImageFileReader< itkImage_t > itkImageReader_t;


    itkImageList_t imageSeries;

    int n = imageNameList.size();
    for (int i = 0; i < n; ++i)
      {
        std::string thisName = imageNameList[i];

        typename itkImageReader_t::Pointer reader = itkImageReader_t::New();
        reader->SetFileName(thisName);

        itkImagePointer_t img;

        try
          {
            reader->Update();
            img = reader->GetOutput();
          }
        catch ( itk::ExceptionObject &err)
          {
            std::cerr<< "ExceptionObject caught !" << std::endl;
            std::cerr<< err << std::endl;
            raise(SIGABRT);
          }


        imageSeries.push_back(img);
      }

    return imageSeries;
  }

  /*
   * readTextLineToListOfString
   */
  template<typename TNull>
  std::vector< std::string > readTextLineToListOfString(const char* textFileName)
  {
    /* The file MUST end with an empty line, then each line will be
       stored as an element in the returned vector object. */


    // Here speed and space is not a big issue, so just let it copy and return stuff.
    std::vector< std::string > listOfStrings;

    std::ifstream f(textFileName);
    std::string thisLine;

    if (f.good())
      {
        while( std::getline(f, thisLine) )
          {
            listOfStrings.push_back(thisLine);
          }
      }
    else
      {
        std::cerr<<"Error: can not open file:"<<textFileName<<std::endl;
        raise(SIGABRT);
      }

    f.close();

    return listOfStrings;
  }


  /**
   * write a component of a vector image
   */
  template< typename itkVectorImage_t >
  void
  writeVectorImage(typename itkVectorImage_t::Pointer img, const char *fileName, int component)
  {
    typedef itk::Image<double, itkVectorImage_t::ImageDimension> itkImage_t;
    typename itkImage_t::Pointer componentImg = itkImage_t::New();
    componentImg->SetRegions(img->GetLargestPossibleRegion() );
    componentImg->Allocate();


    typedef itk::ImageRegionIterator<itkVectorImage_t> VectorIteratorType;
    VectorIteratorType vIter(img, img->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<itkImage_t> IteratorType;
    IteratorType iter(componentImg, componentImg->GetLargestPossibleRegion());

    for (vIter.GoToBegin(), iter.GoToBegin(); !vIter.IsAtEnd(); ++iter, ++vIter)
      {
        iter.Set(vIter.Get()[component]);
      }

    typedef itk::ImageFileWriter< itkImage_t > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(componentImg);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        abort();
      }
  }


  /********************************************************************************
   * write vtkpolydata
   */
  template< typename TNull >
  void writePolydataToXMLFile(const char *fileName, vtkSmartPointer<vtkPolyData> pd)
  {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fileName);
    writer->SetInputData(pd);

//#if VTK_MAJOR_VERSION <= 5
////    writer->SetInput(pd);
//#else
//#endif

    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    writer->SetDataModeToAscii();

    writer->Write();

    return;
  }



  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  readPolyDataFromVtkFile(const char *fileName)
  {
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName( fileName);
    reader->Update();

    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
//    polyData->Update();

    return polyData;
  }

  template< typename TNull >
  vtkSmartPointer<vtkPolyData>
  readPolyDataFromVtpFile(const char *fileName)
  {
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(fileName);
    reader->Update();

    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

    return polyData;
  }


  template< typename TNull >
  std::vector< vtkSmartPointer<vtkPolyData> >
  readPolyDataSeries(const std::vector<std::string>& names)
  {
    std::size_t n = names.size();

    std::vector< vtkSmartPointer<vtkPolyData> > shapeList(n);

    for (std::size_t i = 0; i < n; ++i)
      {
        shapeList[i] = readPolyDataFromVtkFile<char>(names[i].c_str());
      }

    return shapeList;
  }


}// ShapeAnalysis

#endif
