cmake_minimum_required(VERSION 3.5)

project(MFSDA)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MFSDA")
set(EXTENSION_CATEGORY "Stats")
set(EXTENSION_CONTRIBUTORS "John Doe (AnyWare Corp.)")
set(EXTENSION_DESCRIPTION "Slicer Extension of MFSDA_Python")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/MFSDA.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/MFSDA/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(src)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
