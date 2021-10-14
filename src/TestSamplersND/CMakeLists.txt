
include_directories(${PROJECT_SOURCE_DIR}/Screenspace/ZSampler)

add_executable(testIntegrationND testIntegrationND ${PROJECT_SOURCE_DIR}/Screenspace/ZSampler/z.cpp ${PROJECT_SOURCE_DIR}/Screenspace/ZSampler/art2x2-table.cpp ${PROJECT_SOURCE_DIR}/Screenspace/ZSampler/ZTSequence.cpp)
target_link_libraries(testIntegrationND PRIVATE Integrands OpenMP::OpenMP_CXX)

add_executable(testRank testRank.cpp)
target_link_libraries(testRank PRIVATE Integrands OpenMP::OpenMP_CXX)

add_executable(pointsOutput pointsOutput.cpp ${PROJECT_SOURCE_DIR}/Screenspace/ZSampler/z.cpp ${PROJECT_SOURCE_DIR}/Screenspace/ZSampler/art2x2-table.cpp ${PROJECT_SOURCE_DIR}/Screenspace/ZSampler/ZTSequence.cpp)
