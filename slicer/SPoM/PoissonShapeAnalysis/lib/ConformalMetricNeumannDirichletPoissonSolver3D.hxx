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


#ifndef ConformalMetricNeumannDirichletPoissonSolver3D_hxx_
#define ConformalMetricNeumannDirichletPoissonSolver3D_hxx_

#include <vector>
#include <limits>

// itk
#include "itkImage.h"
#include "itkIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_sparse_matrix.h"

// local
#include "ConformalMetricNeumannDirichletPoissonSolver3D.h"
#include "PCGSPDLinearEquationSolver.h"

// dbg
#include "utilitiesIO.h"

namespace ShapeAnalysis
{
  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::update()
  {
    _preprocess();

    _buildIndexMapping();

    _setupEquation();

    _solve();

    m_allDone = true;

    return;
  }


  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_preprocess()
  {
    // TODO add check if images are set
    if (m_inputImage->GetLargestPossibleRegion() != m_DBCMask->GetLargestPossibleRegion())
      {
        std::cerr<<"Error: m_inputImage->GetLargestPossibleRegion() != m_DBCMask->GetLargestPossibleRegion()\n"<<std::flush;
        abort();
      }

    if (m_inputImage->GetLargestPossibleRegion() != m_DBCValueImage->GetLargestPossibleRegion())
      {
        std::cerr<<"Error: m_inputImage->GetLargestPossibleRegion() != m_DBCValueImage->GetLargestPossibleRegion()\n"<<std::flush;
        abort();
      }

    return;
  }

  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_solve()
  {
    /// solve the equation: m_S*m_XUnknown = -m_B*m_XKnown where S is the tri-diagnal matrix(1, 4, 1; ...)
    /// Here, m_S is L_U (in Grady's Random Walker PAMI paper eq8); m_XUnknown is x_U; m_B is B^{\top}; m_XKnown is x_M

    // rhs = -B*m_XKnown
    VNLVectorType rhs(m_numberOfUnknownVoxel);
    m_B.mult(-m_XKnown, rhs);

    // init unknown
    m_XUnknown.set_size(m_numberOfUnknownVoxel);
    m_XUnknown.fill(0.0);

    std::cout<<"Start PCG\n"<<std::flush;

    PCGSPDLinearEquationSolver<FloatType>(m_S, rhs, m_XUnknown, 1e-12);

    m_poissonImage = FloatImageType::New();
    m_poissonImage->SetRegions(m_DBCValueImage->GetLargestPossibleRegion());
    m_poissonImage->Allocate();
    m_poissonImage->FillBuffer(0.0);
    m_poissonImage->CopyInformation(m_DBCValueImage);

    long nx = static_cast<long>(m_DBCValueImage->GetLargestPossibleRegion().GetSize()[0]);
    long ny = static_cast<long>(m_DBCValueImage->GetLargestPossibleRegion().GetSize()[1]);
    long nz = static_cast<long>(m_DBCValueImage->GetLargestPossibleRegion().GetSize()[2]);

    for (long iz = 0; iz < nz; ++iz)
      {
        IndexType idx;
        idx[2] = iz;
        for (long iy = 0; iy < ny; ++iy)
          {
            idx[1] = iy;
            for (long ix = 0; ix < nx; ++ix)
              {
                idx[0] = ix;

                if (m_DBCMask->GetPixel(idx) != 0)
                  {
                    m_poissonImage->SetPixel(idx, m_DBCValueImage->GetPixel(idx));
                  }
                else
                  {
                    m_poissonImage->SetPixel(idx, m_XUnknown[m_indexMap->GetPixel(idx) - m_numberOfDBCVoxel]);
                  }
              }
          }
      }

    return;
  }

  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_buildIndexMapping()
  {
    typedef typename MaskImageType::RegionType RegionType;
    RegionType region = m_DBCMask->GetLargestPossibleRegion();

    m_indexMap = LongImageType::New();
    m_indexMap->SetRegions(region);
    m_indexMap->Allocate();
    m_indexMap->FillBuffer(0);

    //FloatType theMin = std::numeric_limits<FloatType>::min();

    long tmp = 0;

    typedef itk::ImageRegionConstIterator<MaskImageType> ImageRegionConstIterator;
    ImageRegionConstIterator bciter(m_DBCMask, region);

    typedef itk::ImageRegionIterator<LongImageType> ImageRegionIterator;
    ImageRegionIterator iter(m_indexMap, region);

    bciter.GoToBegin();
    iter.GoToBegin();
    for (; !iter.IsAtEnd(); ++iter, ++bciter)
      {
        if (bciter.Get() != 0)
          {
            iter.Set(tmp++);
          }
      }

    m_numberOfDBCVoxel = tmp;

    bciter.GoToBegin();
    iter.GoToBegin();
    for (; !iter.IsAtEnd(); ++iter, ++bciter)
      {
        if (bciter.Get() == 0)
          {
            iter.Set(tmp++);
          }
      }

    m_numberOfUnknownVoxel = tmp - m_numberOfDBCVoxel; // -1 coz tmp has a last increse-by-1 which is not needed.

    // dbg
    std::cout<<"number of all voxels = "<<tmp<<std::endl;

    std::cout<<"number of DBC voxels = "<<m_numberOfDBCVoxel<<",\t number of unknown voxels = "<<m_numberOfUnknownVoxel<<std::endl;

    writeImage<LongImageType>(m_indexMap, "indexMap.nrrd");
    // note: don't check this output using itkSNAP coz it only support short. ImageViewer will work.

    // for (long it = 0; it < m_indexMap.size(); ++it)
    //   {
    //     std::cout<<m_indexMap[it]<<", ";
    //   }
    // std::cout<<std::endl;
    // dbg, end

    return;
  }


  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_setupEquation()
  {
    typedef typename FloatImageType::OffsetType OffsetType;
    OffsetType offsetX;
    offsetX[0] = 1; offsetX[1] = 0; offsetX[2] = 0;

    OffsetType offsetY;
    offsetY[0] = 0; offsetY[1] = 1; offsetY[2] = 0;

    OffsetType offsetZ;
    offsetZ[0] = 0; offsetZ[1] = 0; offsetZ[2] = 1;

    long nx = static_cast<long>(m_DBCMask->GetLargestPossibleRegion().GetSize()[0]);
    long ny = static_cast<long>(m_DBCMask->GetLargestPossibleRegion().GetSize()[1]);
    long nz = static_cast<long>(m_DBCMask->GetLargestPossibleRegion().GetSize()[2]);

    double dx = 1.0;
    double dy = 1.0;
    double dz = 1.0;

    // double dx = static_cast<double>(m_inputImage->GetSpacing()[0]);
    // double dy = static_cast<double>(m_inputImage->GetSpacing()[1]);
    // double dz = static_cast<double>(m_inputImage->GetSpacing()[2]);

    std::cout<<dx<<'\t'<<dy<<'\t'<<dz<<std::endl<<std::flush;

    double dx2resp = 1.0/dx/dx;
    double dy2resp = 1.0/dy/dy;
    double dz2resp = 1.0/dz/dz;

    m_XKnown.set_size(m_numberOfDBCVoxel);
    m_B.set_size(m_numberOfUnknownVoxel, m_numberOfDBCVoxel);

    m_S.set_size(m_numberOfUnknownVoxel, m_numberOfUnknownVoxel);

    for (long iz = 0; iz < nz; ++iz)
      {
        IndexType idx;
        idx[2] = iz;
        for (long iy = 0; iy < ny; ++iy)
          {
            idx[1] = iy;
            for (long ix = 0; ix < nx; ++ix)
              {
                idx[0] = ix;

                long indexInLongVector = m_indexMap->GetPixel(idx);

                if (ix >= 1 && ix < nx-1 && iy >= 1 && iy < ny-1 && iz >= 1 && iz < nz-1)
                  {
                    /// not on the image border

                    if (m_DBCMask->GetPixel(idx) != 0)
                      {
                        m_XKnown[indexInLongVector] = m_DBCValueImage->GetPixel(idx);
                      }
                    else
                      {
                        double degreeTerm = 0.0;

                        {
                          double dfdx = (static_cast<double>(m_inputImage->GetPixel(idx + offsetX)) - static_cast<double>(m_inputImage->GetPixel(idx - offsetX)))/2.0/dx;
                          double conformalFactorX = exp(-m_beta*dfdx*dfdx);
                          double conformalDX2resp = conformalFactorX*dx2resp;

                          //dbg
                          //std::cout<<dfdx<<"                   "<<conformalDX2resp<<"                   ";
                          //dbg, end

                          // x+1
                          if (m_DBCMask->GetPixel(idx + offsetX) != 0)
                            {
                              m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX), conformalDX2resp);
                            }
                          else
                            {
                              m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX) - m_numberOfDBCVoxel, conformalDX2resp);
                            }

                          degreeTerm -= conformalDX2resp;

                          // x-1
                          if (m_DBCMask->GetPixel(idx - offsetX) != 0)
                            {
                              m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX), conformalDX2resp);
                            }
                          else
                            {
                              m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX) - m_numberOfDBCVoxel, conformalDX2resp);
                            }

                          degreeTerm -= conformalDX2resp;
                        }

                        {
                          double dfdy = (static_cast<double>(m_inputImage->GetPixel(idx + offsetY)) - static_cast<double>(m_inputImage->GetPixel(idx - offsetY)))/2.0/dy;
                          double conformalFactorY = exp(-m_beta*dfdy*dfdy);
                          double conformalDY2resp = conformalFactorY*dy2resp;

                          //dbg
                          //std::cout<<dfdy<<"                   "<<conformalDY2resp<<"                   ";
                          //dbg, end

                          // y+1
                          if (m_DBCMask->GetPixel(idx + offsetY) != 0)
                            {
                              m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY), conformalDY2resp);
                            }
                          else
                            {
                              m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY) - m_numberOfDBCVoxel, conformalDY2resp);
                            }

                          degreeTerm -= conformalDY2resp;

                          // y-1
                          if (m_DBCMask->GetPixel(idx - offsetY) != 0)
                            {
                              m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY), conformalDY2resp);
                            }
                          else
                            {
                              m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY) - m_numberOfDBCVoxel, conformalDY2resp);
                            }

                          degreeTerm -= conformalDY2resp;
                        }

                        {
                          double dfdz = (static_cast<double>(m_inputImage->GetPixel(idx + offsetZ)) - static_cast<double>(m_inputImage->GetPixel(idx - offsetZ)))/2.0/dz;
                          double conformalFactorZ = exp(-m_beta*dfdz*dfdz);
                          double conformalDZ2resp = conformalFactorZ*dz2resp;

                          //dbg
                          //std::cout<<dfdz<<"                   "<<conformalDZ2resp<<'\n';
                          //dbg, end


                          // z+1
                          if (m_DBCMask->GetPixel(idx + offsetZ) != 0)
                            {
                              m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ), conformalDZ2resp);
                            }
                          else
                            {
                              m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ) - m_numberOfDBCVoxel, conformalDZ2resp);
                            }
                          degreeTerm -= conformalDZ2resp;

                          // z-1
                          if (m_DBCMask->GetPixel(idx - offsetZ) != 0)
                            {
                              m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ), conformalDZ2resp);
                            }
                          else
                            {
                              m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ) - m_numberOfDBCVoxel, conformalDZ2resp);
                            }
                          degreeTerm -= conformalDZ2resp;
                        }

                        m_S.put(indexInLongVector - m_numberOfDBCVoxel, indexInLongVector - m_numberOfDBCVoxel, degreeTerm);
                      }
                  }
                else
                  {
                    /// On the image border

                    if (m_DBCMask->GetPixel(idx) != 0)
                      {
                        /// come in here only when DBC contains some image border region
                        m_XKnown[indexInLongVector] = m_DBCValueImage->GetPixel(idx);
                      }
                    else
                      {
                        double degreeOfThisPoint = 0.0;

                        // x+1
                        if (ix + 1 < nx)
                          {
                            double dfdx = (static_cast<double>(m_inputImage->GetPixel(idx + offsetX)) - static_cast<double>(m_inputImage->GetPixel(idx)))/dx;
                            double conformalFactorX = exp(-m_beta*dfdx*dfdx);
                            double conformalDX2resp = conformalFactorX/dx;

                            if (m_DBCMask->GetPixel(idx + offsetX) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX), conformalDX2resp);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX) - m_numberOfDBCVoxel, conformalDX2resp);
                              }

                            degreeOfThisPoint -= conformalDX2resp;
                          }


                        // x-1
                        if (ix - 1 >= 0)
                          {
                            double dfdx = (static_cast<double>(m_inputImage->GetPixel(idx)) - static_cast<double>(m_inputImage->GetPixel(idx - offsetX)))/dx;
                            double conformalFactorX = exp(-m_beta*dfdx*dfdx);
                            double conformalDX2resp = conformalFactorX/dx;

                            if (m_DBCMask->GetPixel(idx - offsetX) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX), conformalDX2resp);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX) - m_numberOfDBCVoxel, conformalDX2resp);
                              }

                            degreeOfThisPoint -= conformalDX2resp;
                          }

                        // y+1
                        if (iy + 1 < ny)
                          {
                            double dfdy = (static_cast<double>(m_inputImage->GetPixel(idx + offsetY)) - static_cast<double>(m_inputImage->GetPixel(idx)))/dy;
                            double conformalFactorY = exp(-m_beta*dfdy*dfdy);
                            double conformalDY2resp = conformalFactorY/dy;

                            if (m_DBCMask->GetPixel(idx + offsetY) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY), conformalDY2resp);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY) - m_numberOfDBCVoxel, conformalDY2resp);
                              }
                            degreeOfThisPoint -= conformalDY2resp;
                          }


                        // y-1
                        if (iy - 1 >= 0)
                          {
                            double dfdy = (static_cast<double>(m_inputImage->GetPixel(idx)) - static_cast<double>(m_inputImage->GetPixel(idx - offsetY)))/dy;
                            double conformalFactorY = exp(-m_beta*dfdy*dfdy);
                            double conformalDY2resp = conformalFactorY/dy;

                            if (m_DBCMask->GetPixel(idx - offsetY) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY), conformalDY2resp);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY) - m_numberOfDBCVoxel, conformalDY2resp);
                              }
                            degreeOfThisPoint -= conformalDY2resp;
                          }

                        // z+1
                        if (iz + 1 < nz)
                          {
                            double dfdz = (static_cast<double>(m_inputImage->GetPixel(idx + offsetZ)) - static_cast<double>(m_inputImage->GetPixel(idx)))/dz;
                            double conformalFactorZ = exp(-m_beta*dfdz*dfdz);
                            double conformalDZ2resp = conformalFactorZ/dz;

                            if (m_DBCMask->GetPixel(idx + offsetZ) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ), conformalDZ2resp);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ) - m_numberOfDBCVoxel, conformalDZ2resp);
                              }
                            degreeOfThisPoint -= conformalDZ2resp;
                          }

                        // z-1
                        if (iz - 1 >= 0)
                          {
                            double dfdz = (static_cast<double>(m_inputImage->GetPixel(idx)) - static_cast<double>(m_inputImage->GetPixel(idx - offsetZ)))/dz;
                            double conformalFactorZ = exp(-m_beta*dfdz*dfdz);
                            double conformalDZ2resp = conformalFactorZ/dz;

                            if (m_DBCMask->GetPixel(idx - offsetZ) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ), conformalDZ2resp);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ) - m_numberOfDBCVoxel, conformalDZ2resp);
                              }
                            degreeOfThisPoint -= conformalDZ2resp;
                          }

                        m_S.put(indexInLongVector - m_numberOfDBCVoxel, indexInLongVector - m_numberOfDBCVoxel, degreeOfThisPoint);
                      }
                  }
              }
          }
      }


    // //dbg
    // for (long ir = 0; ir < m_numberOfUnknownVoxel; ++ir)
    //   {
    //     for (long ic = 0; ic < m_numberOfDBCVoxel; ++ic)
    //       {
    //         std::cout<<m_B(ir, ic)<<" ";
    //       }
    //     std::cout<<std::endl;
    //   }
    // //dbg, end

    return;
  }

  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  typename ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::FloatImagePointer
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::getPoissonImage()
  {
    if (m_allDone)
      {
        return m_poissonImage;
      }
    else
      {
        std::cerr<<"Error: not done.\n";
        abort();
      }
  }


  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::setBoundaryConditionMask(MaskImagePointer bcMask)
  {
    m_DBCMask = bcMask;
  }

  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::setBoundaryConditionValue(BCValueImagePointer bcValue)
  {
    m_DBCValueImage = bcValue;
  }


  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::setInputImage(const InputImageType* inputImage)
  {
    m_inputImage = inputImage;

    return;
  }

  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::setBeta(double b)
  {
    if (b > 0)
      {
        m_beta = b;
      }
    else
      {
        std::cout<<"Error: beta must be > 0. Use default 1.0\n"<<std::flush;
      }

    return;
  }



  template< typename TInputImageType, typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  ConformalMetricNeumannDirichletPoissonSolver3D<TInputImageType, TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::ConformalMetricNeumannDirichletPoissonSolver3D()
  {
    m_inputImage = 0;

    m_beta = 1.0;
    m_allDone = false;
  }


}// ShapeAnalysis


#endif
