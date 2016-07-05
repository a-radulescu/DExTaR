# boost build

set (BOOST_BOOTSTRAP_CMD ./bootstrap.sh)
set (BOOST_B2_CMD ./b2)
set (TOOLSET gcc)

# set libraries to build
set (BOOST_LIBS_OPTION --with-system --with-filesystem --with-regex)

ExternalProject_Add(boost
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/boost_1_54_0/
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ${BOOST_BOOTSTRAP_CMD} --prefix=${CMAKE_SOURCE_DIR}/boost_1_54_0/
  BUILD_COMMAND ${BOOST_B2_CMD} install --toolset=${TOOLSET} ${BOOST_LIBS_OPTION}
  INSTALL_COMMAND ""
)

# see also : https://gist.github.com/FlorianWolters/11225791