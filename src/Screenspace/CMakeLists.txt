

###########################################################
################### SCREENSPACE STUFF #####################
## ScreenSpace Laurent's Xcode
add_executable(testSamplerScreenspaceLaurent testSamplerScreenspaceLaurent.cpp)
include_directories(${PROJECT_SOURCE_DIR}/samplerCPP-sig19talk)
add_test(testSamplerScreenspaceLaurent testSamplerScreenspaceLaurent)

add_executable(generateScreenspaceLaurent generateScreenspaceLaurent.cpp)
add_test(generateScreenspaceLaurent generateScreenspaceLaurent)

add_executable(demoOwenPlus2D demoOwenPlus2D.cpp)
target_link_libraries(demoOwenPlus2D PRIVATE OpenMP::OpenMP_CXX)

add_executable(generateScreenspaceSobol++ generateScreenspaceSobol++.cpp)
target_link_libraries(generateScreenspaceSobol++ PRIVATE OpenMP::OpenMP_CXX Integrands)

## tests for the screenspace
add_subdirectory(tests_screenspace)
