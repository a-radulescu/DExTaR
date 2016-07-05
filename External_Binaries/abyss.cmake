# abyss build
set(abyss_PREFIX ${CMAKE_SOURCE_DIR}/abyss)
set(abyss_INSTALL_DIR ${CMAKE_SOURCE_DIR}/)
set(sparsehash_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../External_Libraries/sparsehash/src/sparsehash/include)

set(ABYSS_AUTOGEN_CMD "./autogen.sh")
set(ABYSS_CONFIGURE_CMD "./configure")

message("abyss_PREFIX='${abyss_PREFIX}'")
message("abyss_INSTALL_DIR='${abyss_INSTALL_DIR}'")

set(LIB_DIR ${CMAKE_SOURCE_DIR}/../External_Libraries/)
set(BOOST_ROOT ${LIB_DIR}/boost_1_54_0)
find_package(Boost 1.54 COMPONENTS system filesystem regex REQUIRED)

ExternalProject_Add(abyss
  PREFIX ${abyss_PREFIX}
  SOURCE_DIR ${abyss_PREFIX}
  CONFIGURE_COMMAND cd ${abyss_PREFIX} && ${abyss_PREFIX}/autogen.sh && ${abyss_PREFIX}/${ABYSS_CONFIGURE_CMD} --with-mpi --with-boost=${BOOST_ROOT} --prefix=${abyss_PREFIX} CPPFLAGS=-I${sparsehash_INCLUDE_DIR}

#./configure --with-mpi --with-boost=/home/vintache-d/Logiciels/COMBI/DExTaR/External_Libraries/boost_1_54_0 --prefix=/home/vintache-d/Logiciels/COMBI/DExTaR/External_Binaries/abyss CPPFLAGS=-I/home/vintache-d/Logiciels/COMBI/DExTaR/External_Libraries/sparsehash/src/sparsehash/include


  BUILD_COMMAND cd ${abyss_PREFIX} && make
  INSTALL_COMMAND ""
)

target_link_libraries(abyss ${Boost_LIBRARIES})
