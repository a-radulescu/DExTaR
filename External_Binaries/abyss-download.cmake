# abyss build

execute_process(
  COMMAND git clone --branch 1.3.7 "https://github.com/bcgsc/abyss/"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

#  CONFIGURE_COMMAND ${abyss_PREFIX}/src/abyss/${ABYSS_CONFIGURE_CMD} --prefix=${abyss_PREFIX}/src/abyss
#  INSTALL_COMMAND ""
#)


