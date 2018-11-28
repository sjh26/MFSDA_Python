cmake_minimum_required(VERSION 2.8.9)

project(MFSDA)


set(MFSDA_VERSION_MAJOR 1)
set(MFSDA_VERSION_MINOR 0)
set(MFSDA_VERSION_PATCH 1)


option(MFSDA_BUILD_SLICER_EXTENSION OFF)

if(MFSDA_BUILD_SLICER_EXTENSION)
  #-----------------------------------------------------------------------------
  # Extension meta-information
  set(EXTENSION_HOMEPAGE "")
  set(EXTENSION_CATEGORY "Statistics")
  set(EXTENSION_CONTRIBUTORS "Mateo Lopez (University of North Carolina), Juan Carlos Prieto (University of North Carolina) ")
  set(EXTENSION_DESCRIPTION " ")
  set(EXTENSION_ICONURL "")
  set(EXTENSION_SCREENSHOTURLS "")
  set(EXTENSION_DEPENDS )

  #-----------------------------------------------------------------------------
  # Extension dependencies
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})

  add_subdirectory(MFSDA)
else()
  find_package(VTK REQUIRED)
  find_package(SlicerExecutionModel REQUIRED)
endif()

#-----------------------------------------------------------------------------
# Extension modules

#add_subdirectory(CLI/ComputeMeanShapes)
#add_subdirectory(CLI/linearSubdivision)
## NEXT_MODULE


#-----------------------------------------------------------------------------
if(MFSDA_BUILD_SLICER_EXTENSION)
  include(${Slicer_EXTENSION_CPACK})
endif()
