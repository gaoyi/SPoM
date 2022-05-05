cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SPoM)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/gaoyi/SPoM")
set(EXTENSION_CATEGORY "Shape analysis")
set(EXTENSION_CONTRIBUTORS "Yi Gao (Shenzhen University), Sylvain Bouix (BWH)")
set(EXTENSION_DESCRIPTION "This extension uses the Signed Poisson equations to perform statistical shape analysis.")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/SPoM.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/SPoM/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PoissonShapeAnalysis)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
