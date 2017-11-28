cmake_minimum_required(VERSION 3.5)

#-----------------------------------------------------------------------------
set(MODULE_NAME MFSDA)

option(MFSDA_BUILD_SLICER_EXTENSION OFF)
if( MFSDA_BUILD_SLICER_EXTENSION )

  set(EXTENSION_HOMEPAGE "")
  set(EXTENSION_CATEGORY "Multivariate functional Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "")
  set(EXTENSION_DESCRIPTION "")
  set(EXTENSION_ICONURL "")
  set(EXTENSION_SCREENSHOTURLS "")
  set(EXTENSION_DEPENDS "NA") # Specified as a space separated list or 'NA' if any
  set(EXTENSION_BUILD_SUBDIRECTORY .)
  #unsetForSlicer( NAMES QT_QMAKE_EXECUTABLE SlicerExecutionModel_DIR ITK_DIR VTK_DIR CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS CMAKE_C_FLAGS ITK_LIBRARIES )
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})

  #-----------------------------------------------------------------------------
  set(MODULE_PYTHON_SCRIPTS
    ${MODULE_NAME}.py
    )

  set(MODULE_PYTHON_RESOURCES
    Resources/Icons/${MODULE_NAME}.png
    )

  #-----------------------------------------------------------------------------
  slicerMacroBuildScriptedModule(
    NAME ${MODULE_NAME}
    SCRIPTS ${MODULE_PYTHON_SCRIPTS}
    RESOURCES ${MODULE_PYTHON_RESOURCES}
    WITH_GENERIC_TESTS
    )

  #-----------------------------------------------------------------------------

endif()

if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  #add_subdirectory(Testing)
endif()
