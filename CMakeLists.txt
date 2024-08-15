cmake_minimum_required (VERSION 3.22)

#  To consume the panman library, use:
#
#   FetchContent_Declare(panman
#	   GIT_REPOSITORY https://github.com/TurakhiaLab/pagenome-mat
#	   GIT_TAG [tag]
#	   FIND_PACKAGE_ARGS [...] # to try to find local install first	
#	)
#
#   FetchContent_MakeAvailable(panman)
#
#   target_link_libraries(my_target PRIVATE panman)

include(FetchContent)

project(panman
	VERSION 0.1.0
	DESCRIPTION "Pangenome mutation-annotated trees"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17) # could bump to 20?

set(CMAKE_IGNORE_PREFIX_PATH $ENV{CONDA_PREFIX}/lib/libz.so.1 /usr/lib/tbb)



include(cmake/dependencies.cmake)

FetchContent_MakeAvailable(Boost TBB spoa protobuf jsoncpp)

# pan-mat library
add_library(panman STATIC src/panmanUtils.cpp)

# application binary 
add_executable(panman-utils src/panmanUtils.cpp)


# includes
target_include_directories(panman PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# populate dependencies

# run protoc compiler
include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
set(panman_SRCS panman.pb.cc)
set(panman_PROTOS panman.proto)
target_sources(panman PRIVATE ${panman_SRCS} ${panman_PROTOS})

protobuf_generate(TARGET panman)

# linking
target_link_libraries(panman PUBLIC
	Boost::program_options Boost::filesystem Boost::iostreams Boost::system
	protobuf::libprotobuf
	${TBB_IMPORTED_TARGETS}
	$<TARGET_NAME_IF_EXISTS:jsoncpp_static> $<TARGET_NAME_IF_EXISTS:jsoncpp_lib>
	spoa
)
target_link_libraries(panman-utils PRIVATE panman)


# exporting:
#   - `make` will create build/panman-config.cmake (if make is run)
#   - {install path}/cmake/panman-config.cmake (if make install is run)

set_target_properties(
  panman
  PROPERTIES
  EXPORT_NAME panman-targets
  OUTPUT_NAME panman  # optional: makes the file libTarget.so on disk
)

install(TARGETS panman-utils EXPORT panman-targets)

install(EXPORT panman-targets
    FILE panman-targets.cmake
    NAMESPACE panman::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/panman
)

include(CMakePackageConfigHelpers)

include(GNUInstallDirs)

set(INCLUDE_INSTALL_DIR 
    CACHE PATH "Location of header files" )

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/panmanConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/panman-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/panman
	PATH_VARS INCLUDE_INSTALL_DIR
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/panman-config.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/panman"
)
