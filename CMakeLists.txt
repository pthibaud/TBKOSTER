# For MS windows, try this cmake -G "MinGW Makefiles" -DCMAKE_SH="CMAKE_SH-NOTFOUND" ..
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif()

#cmake_policy(SET CMP0011 NEW)
#cmake_policy(SET CMP0012 NEW)
#cmake_policy(SET CMP0048 NEW)
#cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0074 NEW)

#if(UNIX)
#  if(APPLE OR MINGW)
#    cmake_minimum_required(VERSION 3.9)
#  elseif(NOT APPLE)
#    cmake_minimum_required(VERSION 3.5)
#  endif()
#endif()

cmake_minimum_required(VERSION 3.13)

set(PROJECT_NAME TBKOSTER)
project(${PROJECT_NAME} VERSION 1.9.4.0 DESCRIPTION "Tight-Binding Kernel of Objects for SpinTronics Enhanced Research")
message(STATUS "Preparing ${PROJECT_NAME}-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}.${PROJECT_VERSION_TWEAK} ...")
message(STATUS "Running on ${CMAKE_SYSTEM_NAME} version ${CMAKE_SYSTEM_VERSION}")

set(TBKOSTER_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}.${PROJECT_VERSION_TWEAK}")

# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of TBKOSTER: e.g. "cmake
# -DDEFAULT=0")
set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

opt(DOXYGEN "Enable Doxygen engine for source documentation" OFF)
opt(FLUSH "Enable manually flush files" ON)
opt(LATEX "Enable LaTeX engine for documentation" ON)
opt(OPENMP "Enable OpenMP support" OFF)
opt(TESTS "Enable tests" OFF)
opt(FORTRAN_STDLIB "Enable fortran-lang stdlib" OFF)

# Find all the .in files
file(GLOB inFiles RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/*.in" "${PROJECT_SOURCE_DIR}/test/*.in")
foreach(infileName ${inFiles})
  string(REGEX REPLACE ".in\$" "" outfileName "${infileName}")
  #message(STATUS "Stamp ${infileName} => ${outfileName}")
  configure_file(${PROJECT_SOURCE_DIR}/${infileName} ${PROJECT_SOURCE_DIR}/${outfileName})
endforeach(infileName)

# add custom target for make to clean the distribution of temporary files
add_custom_target(
    distclean
    COMMAND ${CMAKE_COMMAND} -E remove -f "${PROJECT_SOURCE_DIR}/src/*.f90"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${PROJECT_SOURCE_DIR}/test/*.f90"
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

# Enable Fortran language
enable_language(Fortran)
# Force the use of the local version of the compilers
#find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} gfortran PATHS ENV PATH NO_DEFAULT_PATH)
message(STATUS "Using Fortran compiler ${CMAKE_Fortran_COMPILER_ID} Fortran ${CMAKE_Fortran_COMPILER_VERSION}")

# Fortran module and output directories
set(fortran_module_directory "${CMAKE_BINARY_DIR}/modules" CACHE PATH "directory for fortran module files")
file(MAKE_DIRECTORY ${fortran_module_directory})
list(APPEND CMAKE_Fortran_MODULE_DIRECTORY ${fortran_module_directory})
set(fortran_output_directory "${CMAKE_BINARY_DIR}/bin" CACHE PATH "directory for executable file")
file(MAKE_DIRECTORY ${fortran_output_directory})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${fortran_output_directory})

# Enable the preprocessor globally
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
# Generate position independent code (PIC) suitable for use in shared libraries
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")  
set(CMAKE_VERBOSE_MAKEFILE off)

if(ENABLE_FORTRAN_STDLIB)
  # fortran-lang stdlib https://github.com/fortran-lang/stdlib
  find_package("fortran_stdlib" REQUIRED)
  list(APPEND CONFIG_OPTIONS "FORTRAN_STDLIB")
endif(ENABLE_FORTRAN_STDLIB)

# add the directory of the sources
add_subdirectory(src)

# add the directory of the tools
add_subdirectory(tools)

if(ENABLE_FLUSH)
  list(APPEND CONFIG_OPTIONS "FLUSH")
endif(ENABLE_FLUSH)

if(ENABLE_TESTS)
  # unitary tests
  add_subdirectory(unit_tests)
  list(APPEND CONFIG_OPTIONS "TESTS")
endif(ENABLE_TESTS)

if(ENABLE_LATEX)
  add_subdirectory(doc)
endif(ENABLE_LATEX)

if(ENABLE_DOXYGEN)
  # Doxygen documentation
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    message("-- Found DOXYGEN")
    # Add a target to generate API documentation with Doxygen
    set(DOXYGEN_INPUT ${CMAKE_SOURCE_DIR}/doc/Doxyfile)
    add_custom_target(doc
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM)
  else()
    message("Doxygen needs to be installed to generate the Doxygen documentation")
  endif()
  list(APPEND CONFIG_OPTIONS "DOXYGEN")
endif(ENABLE_DOXYGEN)

list(LENGTH CONFIG_OPTIONS SIZE)
if(SIZE GREATER 0)
	list(SORT CONFIG_OPTIONS)
else()
	list(APPEND CONFIG_OPTIONS "NONE")
endif()

foreach(OPT ${CONFIG_OPTIONS})
	set(TBKOSTER_CONFIG_OPTIONS "${TBKOSTER_CONFIG_OPTIONS} ${OPT}")
endforeach()

message(STATUS "")
message(STATUS "TBKOSTER ${TBKOSTER_VERSION} has been configured for ${CMAKE_HOST_SYSTEM}")
message(STATUS "")
message(STATUS " * Build options: " ${TBKOSTER_CONFIG_OPTIONS})
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * Fortran compiler: " ${CMAKE_Fortran_COMPILER} " Release: " ${CMAKE_Fortran_COMPILER_VERSION})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
