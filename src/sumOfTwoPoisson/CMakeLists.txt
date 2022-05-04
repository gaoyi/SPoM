###############################################################################
#                                                                             #
#  Copyright (c) 2010, Yi Gao                                                 #
#  gaoyi@gatech.edu                                                           #
#                                                                             #
#                                                                             #
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
#  DEALINGS IN THE SOFTWARE.                                                  #
#                                                                             #
#  See the README.md and COPYING files for usage and copyright information.   #
#                                                                             #
###############################################################################


#include_directories( ${PROJECT_SOURCE_DIR}/sumOfTwoPoisson/lib )
include_directories( lib )

################################################################################
# exe's
option(build_mainSumOfTwoPoissonShapeAnalysis "build mainSumOfTwoPoissonShapeAnalysis?" ON)
if (build_mainSumOfTwoPoissonShapeAnalysis)
    set(cexx mainSumOfTwoPoissonShapeAnalysis)
    # add_executable(${cexx} ${cexx}.cxx )
    # target_link_libraries(${cexx} alglib localLibrary ${namicLibraries})

    add_executable(${cexx} ${cexx}.cxx ${headerAndLibDir}/itkVtkMeshConversion.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_mainSumOfTwoPoissonShapeAnalysis)



option(build_testShapeSumOfTwoPoissonFilter "build testShapeSumOfTwoPoissonFilter?" ON)
if (build_testShapeSumOfTwoPoissonFilter)
    set(cexx testShapeSumOfTwoPoissonFilter)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_testShapeSumOfTwoPoissonFilter)



option(build_testConformalMetricNeumannDirichletPoissonSolver3D "build testConformalMetricNeumannDirichletPoissonSolver3D?" ON)
if (build_testConformalMetricNeumannDirichletPoissonSolver3D)
    set(cexx testConformalMetricNeumannDirichletPoissonSolver3D)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_testConformalMetricNeumannDirichletPoissonSolver3D)




option(build_mainFlipVTKPolydata "build mainFlipVTKPolydata?" ON)
if (build_mainFlipVTKPolydata)
    set(cexx mainFlipVTKPolydata)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_mainFlipVTKPolydata)




option(build_mainProbeVolumeWithMesh "build mainProbeVolumeWithMesh?" ON)
if (build_mainProbeVolumeWithMesh)
    set(cexx mainProbeVolumeWithMesh)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_mainProbeVolumeWithMesh)



option(build_mainRegionCoverRate "build mainRegionCoverRate?" ON)
if (build_mainRegionCoverRate)
    set(cexx mainRegionCoverRate)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_mainRegionCoverRate)




option(build_mainProbeVolumeWithMeshAndComputeDice "build mainProbeVolumeWithMeshAndComputeDice?" ON)
if (build_mainProbeVolumeWithMeshAndComputeDice)
    set(cexx mainProbeVolumeWithMeshAndComputeDice)
    add_executable(${cexx} ${cexx}.cxx)
    target_link_libraries(${cexx} alglib ${namicLibraries})
endif (build_mainProbeVolumeWithMeshAndComputeDice)
