project(cmake-openmp)
cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add extra include directories
include_directories(/usr/local/include)


# add extra lib directories
link_directories(/usr/local/opt/libomp/lib /usr/local/lib)

if(APPLE)
    message(">>>>>>>>>>>>>>>>>>>>>>>>>>  Setting OpenMP flags on MacOs. Assuming `brew install libomp`")
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

#SET(CMAKE_CXX_FLAGS "-lm")

set(Integrands ../Integration/Cont2D.cpp ../Integration/Gauss1D.cpp ../Integration/Gauss2D.cpp ../Integration/Gauss3D.cpp ../Integration/Gauss4D.cpp ../Integration/Gauss5D.cpp ../Integration/Gauss6D.cpp ../Integration/Gauss8D.cpp ../Integration/Gauss10D.cpp ../Integration/Gauss12D.cpp ../Integration/GaussIso2D.cpp ../Integration/HeaviBell2D.cpp ../Integration/HeaviCont2D.cpp ../Integration/HeaviGauss2D.cpp ../Integration/Heaviside2D.cpp ../Integration/Heaviside3D.cpp ../Integration/Heaviside4D.cpp ../Integration/Heaviside5D.cpp ../Integration/Heaviside6D.cpp ../Integration/Smooth2D.cpp ../Integration/Smooth3D.cpp ../Integration/Smooth4D.cpp ../Integration/SmoothBump2D.cpp)

find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g -DLOGLEVEL=3 -DPFAEDLE_DBG=1")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS} -DLOGLEVEL=2")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS} -DLOGLEVEL=2")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -g -DLOGLEVEL=3")

#######################################################################################

add_executable(owen owen.cpp)
target_link_libraries(owen PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS owen DESTINATION ~/bin)

add_executable(getUniformNDpts getUniformNDpts.cpp)
target_link_libraries(getUniformNDpts PRIVATE OpenMP::OpenMP_CXX)
install(TARGETS getUniformNDpts DESTINATION ~/bin)

#add_executable(owen64 owen64.cpp)
#target_link_libraries(owen64 PRIVATE OpenMP::OpenMP_CXX)
#install(TARGETS owen64 DESTINATION ~/bin)

#add_executable(owen64global owen64global.cpp)
#target_link_libraries(owen64global PRIVATE OpenMP::OpenMP_CXX)
#install(TARGETS owen64global DESTINATION ~/bin)

SET(CMAKE_INSTALL_PREFIX $ENV{HOME})
