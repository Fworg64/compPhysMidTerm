cmake_minimum_required (VERSION 2.8.1)
project(PROGRAM)
set(CMAKE_BUILD_TYPE Debug)

set(GRAPHLIB_HOME "/home/fworg64/lib")
#set(GSL_LIBS "/usr/local/lib/gsl/lib/libgsl.so /usr/local/lib/gsl/lib/libgslcblas.so m")

include_directories("~/include")
include_directories("/usr/local/lib/gsl/include")

add_executable(programexec src/program.c)

target_link_libraries(programexec 
  ${GRAPHLIB_HOME}/libgraph.a 
  m 
  X11 
  "/usr/local/lib/gsl/lib/libgsl.so" 
  "/usr/local/lib/gsl/lib/libgslcblas.so"
)


add_executable(gsltest src/gsl_example.c)
#link_directories("/usr/local/lib/gsl/lib")
target_link_libraries(gsltest 
  "/usr/local/lib/gsl/lib/libgsl.so" 
  "/usr/local/lib/gsl/lib/libgslcblas.so"
  m
)