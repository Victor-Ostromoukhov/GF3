set(INTEGRANDS ../../Integration/Gauss2D.cpp)

if(WIN32)
  set(INTEGRANDS ${INTEGRANDS} ../../Tools/drand48.cpp)
endif()

add_executable(integrationErrorPerPixel integrationErrorPerPixel.cpp ${INTEGRANDS})
