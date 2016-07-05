cmake_minimum_required (VERSION 2.8)
project (External_Libraries)

set(LIB_DIR "${PROJECT_SOURCE_DIR}/External_Libraries")

include(ExternalProject)
include(boost-download.cmake)
include(boost.cmake)
include(sparsehash.cmake)
