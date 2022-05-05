###############################################################################
#                                                                             #
#  This program is distributed in the hope that it will be useful, but        #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY #
#  or FITNESS FOR A PARTICULAR PURPOSE.                                       #
#                                                                             #
#  Please contact the author for reusing or redistributing this program.      #
#                                                                             #
#                                                  Copyright (c) 2013, Yi Gao #
#                                                            gaoyi@gatech.edu #
#                                                                             #
###############################################################################


cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(projName "ShapeAnalysis")
project(${projName})
# INCLUDE_REGULAR_EXPRESSION("^.*$")


set(CMAKE_CXX_FLAGS "-fopenmp")
set(CMAKE_C_FLAGS "-fopenmp")


################################################################################
# VTK
find_package(VTK COMPONENTS REQUIRED
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  QUIET
            )

if (NOT VTK_FOUND)
  message("Skipping Dodecahedron: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")


################################################################################
################################
# include ITK, to avoid setting it everytime, do:
# set the environment variable ITK_DIR to the dir containing ITKConfig.cmake in .bashrc
#set(ITK_DIR $ENV{ITK_DIR})

#find_package (ITK REQUIRED ITKCommon ITKVtkGlue)  # this does not work
find_package (ITK REQUIRED)
# This askes that when building ITK, should turn itkVTKGlue on. The source of them is at InsightToolkit-4.1.0/Modules/Bridge/VtkGlue/include. There are in fact more files in there than in the separate itkVTKGlue package downloaded from ITK/example site.
if (ITK_FOUND)
   include(${ITK_USE_FILE})
endif(ITK_FOUND)


# if (VTK_FOUND)
#    include(${VTK_USE_FILE})
# endif(VTK_FOUND)

################################################################################
# find gsl path
find_path(GSL_INSTALLATION_PATH include/gsl/gsl_math.h REQUIRED)
if(NOT GSL_INSTALLATION_PATH)
   message(FATAL_ERROR "Please enter the installation path of gsl, i.e., the one after --prefix= when configuring gsl")
endif(NOT GSL_INSTALLATION_PATH)

include_directories(${GSL_INSTALLATION_PATH}/include)

set(gslLibs ${GSL_INSTALLATION_PATH}/lib/libgsl.a)

set(itkLibraries ITKQuadEdgeMesh ${ITK_LIBRARIES})

set(namicLibraries ${VTK_LIBRARIES} ${itkLibraries})


################################################################################
# alglib
set(ALGLIB_PATH ${PROJECT_SOURCE_DIR}/lib/alglib/src)

find_path(ALGLIB_PATH ap.h)
if(NOT ALGLIB_PATH)
   message(FATAL_ERROR "Please enter the path of alglib")
endif(NOT ALGLIB_PATH)

include_directories(${ALGLIB_PATH})
#lind_directories(../)

add_library(alglib
  ${ALGLIB_PATH}/alglibinternal.cpp
  ${ALGLIB_PATH}/alglibmisc.cpp
  ${ALGLIB_PATH}/ap.cpp
  ${ALGLIB_PATH}/dataanalysis.cpp
  ${ALGLIB_PATH}/diffequations.cpp
  ${ALGLIB_PATH}/fasttransforms.cpp
  ${ALGLIB_PATH}/integration.cpp
  ${ALGLIB_PATH}/interpolation.cpp
  ${ALGLIB_PATH}/linalg.cpp
  ${ALGLIB_PATH}/optimization.cpp
  ${ALGLIB_PATH}/solvers.cpp
  ${ALGLIB_PATH}/specialfunctions.cpp
  ${ALGLIB_PATH}/statistics.cpp
  )


################################################################################



################################################################################
# local headers and libs
set(headerAndLibDir ${PROJECT_SOURCE_DIR}/lib)
include_directories(${headerAndLibDir})

add_library(localLibrary ${headerAndLibDir}/itkVtkMeshConversion.cxx)
target_link_libraries(localLibrary alglib ${namicLibraries})


################################################################################
add_subdirectory(sumOfTwoPoisson)
