PROJECT(gkitpath)

cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    add_definitions("-DGK_MACOS")
    add_definitions("-DGL_SILENCE_DEPRECATION")
elseif(WIN32)
  add_definitions("-DWIN32")
  add_definitions("-D_USE_MATH_DEFINES")
  add_definitions("-D_CRT_SECURE_NO_WARNINGS")
  add_definitions("-DNOMINMAX")
  find_package(GLEW REQUIRED)
  include_directories(${GLEW_INCLUDE_DIRS})
elseif(UNIX)
  find_package(GLEW REQUIRED)
  include_directories(${GLEW_INCLUDE_DIRS})
endif()

find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/path")
include_directories("${PROJECT_SOURCE_DIR}/src/gKit")

file(GLOB path_SRC_gKit
    "${PROJECT_SOURCE_DIR}/src/gKit/*.cpp"
)
set(path_SRC
    "path/main.cpp" 
    "path/stbvh_builder.cpp" "path/bvh_builder.cpp" "path/stbvh.cpp" "path/bvh.cpp"
    "path/framebuffer.cpp" "path/texturelib.cpp" "path/materials.cpp" "path/scene.cpp"
    "path/sampler_s19.cpp" "path/sobol_spack.cpp" 
    "../Screenspace/ZSampler/zhash.cpp"
    "../Screenspace/ZSampler/z.cpp" "../Screenspace/ZSampler/zpp.cpp" "../Screenspace/ZSampler/art2x2-table.cpp"
    
)

find_package(SDL2 REQUIRED)
find_package(SDL2_image REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS}  ${SDL2_IMAGE_INCLUDE_DIRS} "/usr/local/include/")

add_executable(path6d  ${path_SRC} ${path_SRC_gKit})
target_link_libraries(path6d PRIVATE OpenMP::OpenMP_CXX OpenGL::GL ${SDL2_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_IMAGE_LIBRARIES})

set(blue_SRC
    "path/blue.cpp" 
    "path/stbvh_builder.cpp" "path/bvh_builder.cpp" "path/stbvh.cpp" "path/bvh.cpp"
    "path/framebuffer.cpp" "path/texturelib.cpp" "path/materials.cpp" "path/scene.cpp"
    "path/sampler_s19.cpp" "path/sobol_spack.cpp" 
)

add_executable(blue6d  ${blue_SRC} ${path_SRC_gKit})
target_link_libraries(blue6d PRIVATE OpenMP::OpenMP_CXX OpenGL::GL ${SDL2_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_IMAGE_LIBRARIES})

add_executable(grid  ${PROJECT_SOURCE_DIR}/src/grid.cpp  ${path_SRC_gKit})
target_link_libraries(grid PRIVATE OpenMP::OpenMP_CXX OpenGL::GL ${SDL2_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_IMAGE_LIBRARIES})

add_executable(stats  ${PROJECT_SOURCE_DIR}/src/stats.cpp  ${path_SRC_gKit})
target_link_libraries(stats PRIVATE OpenMP::OpenMP_CXX OpenGL::GL ${SDL2_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_IMAGE_LIBRARIES})

add_executable(hdrtopng  ${PROJECT_SOURCE_DIR}/src/hdrtopng.cpp  ${path_SRC_gKit})
target_link_libraries(hdrtopng PRIVATE OpenMP::OpenMP_CXX OpenGL::GL ${SDL2_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_IMAGE_LIBRARIES})


add_executable(seeds "src/seeds.cpp")

