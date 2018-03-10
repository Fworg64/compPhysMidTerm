cmake_minimum_required (VERSION 2.8.1)
project(PROGRAM)
set(CMAKE_BUILD_TYPE Debug)

set(GRAPHLIB_HOME "/home/fworg64/lib")
#set(GSL_LIBS "/usr/local/lib/gsl/lib/libgsl.so /usr/local/lib/gsl/lib/libgslcblas.so m")

include_directories("~/include")
include_directories("/usr/local/lib/gsl/include")

include_directories(include)

add_library(renderer src/renderer.c)

target_link_libraries(renderer
  png
  m
  "/usr/local/lib/gsl/lib/libgsl.so" 
  "/usr/local/lib/gsl/lib/libgslcblas.so"
)

############
# Build Test programs
############
add_executable(raygenTest src/tests/cameraRayGenerationTest.c)

target_link_libraries(raygenTest 
  renderer
)
#############
#Build main program
#############
#add_executable(programexec src/main.c)

#target_link_libraries(programexec 
#  ${GRAPHLIB_HOME}/libgraph.a 
#  m 
#  X11 
#  renderer
#)


##############
#example for compiling gsl code
##############
#add_executable(gsltest src/gsl_example.c)
#
#target_link_libraries(gsltest 
#  "/usr/local/lib/gsl/lib/libgsl.so" 
#  "/usr/local/lib/gsl/lib/libgslcblas.so"
#  m
#)
#