
#-----------------------------------------------------------------------------
set(MODULE_NAME PoissonShapeAnalysis)

#-----------------------------------------------------------------------------

#
# SlicerExecutionModel
#
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})


#added by gaoyi ###############################################################################

# include(setupExternalPackages.cmake)


################################################################################
# find gsl path
find_path(GSL_INSTALLATION_PATH include/gsl/gsl_math.h REQUIRED)
if(NOT GSL_INSTALLATION_PATH)
   message(FATAL_ERROR "Please enter the installation path of gsl, i.e., the one after --prefix= when configuring gsl")
endif(NOT GSL_INSTALLATION_PATH)

include_directories(${GSL_INSTALLATION_PATH}/include)

set(gslLibs ${GSL_INSTALLATION_PATH}/lib/libgsl.a)

set(namicLibraries ${VTK_LIBRARIES} ${itkLibraries})


################################################################################
# alglib
#set(ALGLIB_PATH ${PROJECT_SOURCE_DIR}/lib/alglib/src)
set(ALGLIB_PATH lib/alglib/src)

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
# set(headerAndLibDir lib sumOfTwoPoisson)
# include_directories(${headerAndLibDir})
# include_directories(${headerAndLibDir})

#added by gaoyi ###############################################################################


#
# ITK
#
set(${PROJECT_NAME}_ITK_COMPONENTS
    ITKIOImageBase
    ITKCommon
    ITKIOTransformBase
    ITKTransform
    ITKSmoothing
    ITKOptimizers
    ITKQuadEdgeMesh
   )
find_package(ITK 4.6 COMPONENTS ${${PROJECT_NAME}_ITK_COMPONENTS} REQUIRED)
if(ITK_VERSION VERSION_GREATER_EQUAL "5.3")
  foreach(factory_uc IN ITEMS "IMAGEIO" "MESHIO" "TRANSFORMIO")
    set(ITK_NO_${factory_uc}_FACTORY_REGISTER_MANAGER 1)
  endforeach()
else()
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # See Libs/ITKFactoryRegistration/CMakeLists.txt
endif()
include(${ITK_USE_FILE})

set(itkLibraries ITKQuadEdgeMesh ${ITK_LIBRARIES})



#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
    lib
  )


add_library(localLibrary lib/itkVtkMeshConversion.cxx)
target_link_libraries(localLibrary alglib ${ITK_LIBRARIES} ${VTK_LIBRARIES})


set(MODULE_SRCS
#    ${headerAndLibDir}/itkVtkMeshConversion.cxx
  )

set(MODULE_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES} alglib localLibrary ${namicLibraries}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
