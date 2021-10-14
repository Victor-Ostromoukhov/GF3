PROJECT(Sobol++)

cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

option(WITH_GKIT "Build with GKIT/path6d" ON)
option(WITH_GEOGRAM "Build with geogram/OT" OFF)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

if(APPLE)
    message("-- Setting OpenMP flags on MacOs. Assuming `brew install libomp`")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    endif()
endif()
find_package(OpenMP REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/includes)
include_directories(${PROJECT_SOURCE_DIR}/)

add_executable(otDistanceMatrix otDistanceMatrix.cpp)
target_link_libraries(otDistanceMatrix PRIVATE OpenMP::OpenMP_CXX)

set(INTEGRANDS Integration/Cont2D.cpp Integration/Gauss1D.cpp Integration/Gauss2D.cpp Integration/Gauss3D.cpp Integration/Gauss4D.cpp Integration/Gauss5D.cpp Integration/Gauss6D.cpp Integration/Gauss8D.cpp Integration/Gauss10D.cpp Integration/Gauss12D.cpp Integration/GaussIso2D.cpp Integration/HeaviBell2D.cpp Integration/HeaviCont2D.cpp Integration/HeaviGauss2D.cpp Integration/Heaviside2D.cpp Integration/Heaviside3D.cpp Integration/Heaviside4D.cpp Integration/Heaviside5D.cpp Integration/Heaviside6D.cpp Integration/Smooth2D.cpp Integration/Smooth3D.cpp Integration/Smooth4D.cpp Integration/SmoothBump2D.cpp)

if(WIN32)
  set(INTEGRANDS ${INTEGRANDS} Tools/drand48.cpp)
endif()

add_library(Integrands STATIC ${INTEGRANDS} )


set(INTEGRANDS)
if(WIN32)
  set(INTEGRANDS Tools/drand48.cpp)
endif()


include_directories(${PROJECT_SOURCE_DIR}/pcg32-demo)
add_executable(pcg32-demo pcg32-demo/pcg32-demo)
target_link_libraries(pcg32-demo PRIVATE Integrands OpenMP::OpenMP_CXX)

##Geogram
if (WITH_GEOGRAM)
  include(CPM)
  if(WIN32)
    set(VORPALINE_PLATFORM2 Win-vs-dynamic-generic)
  elseif(APPLE)
    set(VORPALINE_PLATFORM2 Darwin-clang-dynamic)
  elseif(UNIX)
    set(VORPALINE_PLATFORM2 Linux64-gcc-dynamic)
  endif()
  set(GEOGRAM_WITH_EXPLORAGRAM ON)
  CPMAddPackage(
    NAME   geogram
    VERSION 1.7.6
    GITHUB_REPOSITORY alicevision/geogram
    OPTIONS
    "GEOGRAM_LIB_ONLY ON"
    "VORPALINE_PLATFORM ${VORPALINE_PLATFORM2}"
    "GEOGRAM_WITH_EXPLORAGRAM ON"
    )

  include_directories(${geogram_SOURCE_DIR}/src/lib/)

  add_executable(otCostCompare otCostCompare.cpp Transport/Wasserstein.cpp Transport/NBallRadonManager.cpp Points/mapping.cpp Samplers/myMath.cpp Tools/my_utility.cpp  Integration/Gauss2D.cpp ${INTEGRANDS})
  target_link_libraries(otCostCompare PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(otCostCompare PRIVATE geogram)
  target_link_directories(otCostCompare PRIVATE ${PROJECT_BINARY_DIR}/lib/)
  target_link_libraries(otCostCompare PRIVATE exploragram)

  find_package(GSL)
  if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    target_link_libraries(otCostCompare PRIVATE ${GSL_LIBRARIES})
  else()
    CPMAddPackage(
      NAME gsl
      GITHUB_REPOSITORY ampl/gsl
      VERSION 2.5.0
      OPTIONS "GSL_DISABLE_TESTS ON"
      )
    target_link_libraries(otCostCompare PRIVATE gsl)
  endif()
endif()

add_executable(viewFunction Integration/viewGaussian.cpp)
target_link_libraries(viewFunction PRIVATE Integrands)

add_executable(testTileGaussian Screenspace/testTileGaussian.cpp)
target_link_libraries(testTileGaussian PRIVATE Integrands OpenMP::OpenMP_CXX)

add_executable(cascadedSobol cascadedSobol.cpp)

## gk2path6D
if (WITH_GKIT)
  add_subdirectory(gk2path6d)
endif()

#Screenspace
add_subdirectory(Screenspace)

#Optimization
add_subdirectory(Optimization)

#testND
add_subdirectory(TestSamplersND)

#vo_tools
add_subdirectory(vo_tools)
