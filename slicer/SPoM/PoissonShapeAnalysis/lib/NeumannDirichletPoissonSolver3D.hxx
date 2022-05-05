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


#ifndef NeumannDirichletPoissonSolver3D_hxx_
#define NeumannDirichletPoissonSolver3D_hxx_

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
#include "NeumannDirichletPoissonSolver3D.h"
#include "PCGSPDLinearEquationSolver.h"

// dbg
#include "utilitiesIO.h"

namespace ShapeAnalysis
{
  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::update()
  {
    _preprocess();

    _buildIndexMapping();

    _setupEquation();

    _solve();

    m_allDone = true;

    return;
  }


  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_preprocess()
  {
    // TODO add check if images are set

    return;
  }

  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_solve()
  {
    // solve the equation: S*m_XUnknown = 1-B*m_XKnown where S is the tri-diagnal matrix(1, 4, 1; ...)

    // rhs = -B*m_XKnown
    VNLVectorType rhs(m_numberOfUnknownVoxel);
    m_B.mult(-m_XKnown, rhs);
    rhs = rhs + 1; // without this 1, it's a Laplace eq

    // init unknown
    m_XUnknown.set_size(m_numberOfUnknownVoxel);
    m_XUnknown.fill(0.0);

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

  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_buildIndexMapping()
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


  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::_setupEquation()
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
                    if (m_DBCMask->GetPixel(idx) != 0)
                      {
                        m_XKnown[indexInLongVector] = m_DBCValueImage->GetPixel(idx);
                      }
                    else
                      {
                        m_S.put(indexInLongVector - m_numberOfDBCVoxel, indexInLongVector - m_numberOfDBCVoxel, -6.0);

                        // x+1
                        if (m_DBCMask->GetPixel(idx + offsetX) != 0)
                          {
                            m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX), 1.0);
                          }
                        else
                          {
                            m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX) - m_numberOfDBCVoxel, 1.0);
                          }

                        // x-1
                        if (m_DBCMask->GetPixel(idx - offsetX) != 0)
                          {
                            m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX), 1.0);
                          }
                        else
                          {
                            m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX) - m_numberOfDBCVoxel, 1.0);
                          }

                        // y+1
                        if (m_DBCMask->GetPixel(idx + offsetY) != 0)
                          {
                            m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY), 1.0);
                          }
                        else
                          {
                            m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY) - m_numberOfDBCVoxel, 1.0);
                          }


                        // y-1
                        if (m_DBCMask->GetPixel(idx - offsetY) != 0)
                          {
                            m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY), 1.0);
                          }
                        else
                          {
                            m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY) - m_numberOfDBCVoxel, 1.0);
                          }

                        // z+1
                        if (m_DBCMask->GetPixel(idx + offsetZ) != 0)
                          {
                            m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ), 1.0);
                          }
                        else
                          {
                            m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ) - m_numberOfDBCVoxel, 1.0);
                          }

                        // z-1
                        if (m_DBCMask->GetPixel(idx - offsetZ) != 0)
                          {
                            m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ), 1.0);
                          }
                        else
                          {
                            m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ) - m_numberOfDBCVoxel, 1.0);
                          }
                      }
                  }
                else
                  {
                    if (m_DBCMask->GetPixel(idx) != 0)
                      {
                        m_XKnown[indexInLongVector] = m_DBCValueImage->GetPixel(idx);
                      }
                    else
                      {
                        double degreeOfThisPoint = 0.0;

                        // x+1
                        if (ix + 1 < nx)
                          {
                            if (m_DBCMask->GetPixel(idx + offsetX) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX), 1.0);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetX) - m_numberOfDBCVoxel, 1.0);
                              }

                            degreeOfThisPoint -= 1.0;
                          }


                        // x-1
                        if (ix - 1 >= 0)
                          {
                            if (m_DBCMask->GetPixel(idx - offsetX) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX), 1.0);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetX) - m_numberOfDBCVoxel, 1.0);
                              }

                            degreeOfThisPoint -= 1.0;
                          }

                        // y+1
                        if (iy + 1 < ny)
                          {
                            if (m_DBCMask->GetPixel(idx + offsetY) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY), 1.0);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetY) - m_numberOfDBCVoxel, 1.0);
                              }
                            degreeOfThisPoint -= 1.0;
                          }


                        // y-1
                        if (iy - 1 >= 0)
                          {
                            if (m_DBCMask->GetPixel(idx - offsetY) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY), 1.0);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetY) - m_numberOfDBCVoxel, 1.0);
                              }
                            degreeOfThisPoint -= 1.0;
                          }

                        // z+1
                        if (iz + 1 < nz)
                          {
                            if (m_DBCMask->GetPixel(idx + offsetZ) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ), 1.0);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx + offsetZ) - m_numberOfDBCVoxel, 1.0);
                              }
                            degreeOfThisPoint -= 1.0;
                          }

                        // z-1
                        if (iz - 1 >= 0)
                          {
                            if (m_DBCMask->GetPixel(idx - offsetZ) != 0)
                              {
                                m_B.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ), 1.0);
                              }
                            else
                              {
                                m_S.put(indexInLongVector - m_numberOfDBCVoxel, m_indexMap->GetPixel(idx - offsetZ) - m_numberOfDBCVoxel, 1.0);
                              }
                            degreeOfThisPoint -= 1.0;
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

  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  typename NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::FloatImagePointer
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::getPoissonImage()
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


  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::setBoundaryConditionMask(MaskImagePointer bcMask)
  {
    m_DBCMask = bcMask;
  }

  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  void
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::setBoundaryConditionValue(BCValueImagePointer bcValue)
  {
    m_DBCValueImage = bcValue;
  }

  template< typename TNeumannDirichletBoundaryConditionMaskType, typename TNeumannDirichletBoundaryConditionValueImageType >
  NeumannDirichletPoissonSolver3D<TNeumannDirichletBoundaryConditionMaskType, TNeumannDirichletBoundaryConditionValueImageType>::NeumannDirichletPoissonSolver3D()
  {
    m_allDone = false;
  }

}// ShapeAnalysis


#endif
