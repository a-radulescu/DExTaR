cmake_minimum_required (VERSION 2.8)
project (DExTaR_External_Binaries)

set(LIB_DIR "${PROJECT_SOURCE_DIR}/External_Binaries")

include(ExternalProject)
include(abyss-download.cmake)
include(abyss.cmake)
