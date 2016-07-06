cmake_minimum_required (VERSION 2.8)
project (DExTaR)

set(CMAKE_CXX_FLAGS "-O3")
set(BOOST_ROOT External_Libraries/boost_1_54_0/)
include_directories(
  External_Libraries/boost_1_54_0/include
  External_Libraries/boost_1_54_0/boost/graph
  )
set(Boost_USE_STATIC_LIBS TRUE)
find_package(Boost 1.54 COMPONENTS system filesystem REQUIRED)

set(DEXTAR_SRC src/DExTaR.cpp)

# define the binary name to compile
add_executable(dextar ${DEXTAR_SRC})
# link to specified library
target_link_libraries(dextar ${Boost_LIBRARIES})
# install binary in specified directory
install(TARGETS dextar RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})