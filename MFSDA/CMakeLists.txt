cmake_minimum_required(VERSION 3.12)

#-----------------------------------------------------------------------------
set(MODULE_NAME MFSDA)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/UI/${MODULE_NAME}.ui
  Resources/Libraries/MFSDA_stat.py
  Resources/Libraries/stat_kernel.py
  Resources/Libraries/stat_lpks.py
  Resources/Libraries/stat_sif.py
  Resources/Libraries/stat_wald_ht.py
  Resources/Libraries/stat_bstrp_pvalue.py
  Resources/Libraries/stat_read_x.py
  Resources/Libraries/stat_lpks_pre_bw.py
  Resources/Libraries/stat_grs.py
  MFSDA_run.py
  MFSDA_run.xml
  MFSDA_createShapes.py
  MFSDA_createShapes.xml
  )




find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  #add_subdirectory(Testing)
endif()
