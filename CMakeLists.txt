cmake_minimum_required (VERSION 3.22)

#  To consume the panmat library, use:
#
#   FetchContent_Declare(panmat
#	   GIT_REPOSITORY https://github.com/TurakhiaLab/pagenome-mat
#	   GIT_TAG [tag]
#	   FIND_PACKAGE_ARGS [...] # to try to find local install first	
#	)
#
#   FetchContent_MakeAvailable(panmat)
#
#   target_link_libraries(my_target PRIVATE panmat)

include(cmake/dependencies.cmake)

project(panmat
	VERSION 0.1.0
	DESCRIPTION "Pangenome mutation-annotated trees"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17) # could bump to 20?

set(CMAKE_IGNORE_PREFIX_PATH $ENV{CONDA_PREFIX}/lib/libz.so.1)

# pan-mat library
add_library(panmat STATIC src/PangenomeMAT.cpp)

# application binary 
add_executable(panmat-utils src/main.cpp src/PangenomeMAT.cpp)

# includes
target_include_directories(panmat PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# populate dependencies
FetchContent_MakeAvailable(Boost TBB spoa protobuf jsoncpp)

# run protoc compiler
include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
set(panmat_SRCS mutation_annotation_test_proto3_optional_new.pb.cc)
set(panmat_PROTOS mutation_annotation_test_proto3_optional_new.proto)
target_sources(panmat PRIVATE ${panmat_SRCS} ${panmat_PROTOS})

protobuf_generate(TARGET panmat)

# linking
target_link_libraries(panmat PUBLIC
	Boost::program_options Boost::filesystem Boost::iostreams Boost::system
	protobuf::libprotobuf
	TBB::tbb
	$<TARGET_NAME_IF_EXISTS:jsoncpp_static> $<TARGET_NAME_IF_EXISTS:jsoncpp_lib>
	spoa
)
target_link_libraries(panmat-utils PRIVATE panmat)


# exporting:
#   - `make` will create build/panmat-config.cmake (if make is run)
#   - {install path}/cmake/panmat-config.cmake (if make install is run)

set_target_properties(
  panmat
  PROPERTIES
  EXPORT_NAME panmat-targets
  OUTPUT_NAME panmat  # optional: makes the file libTarget.so on disk
)

install(TARGETS panmat-utils EXPORT panmat-targets)

install(EXPORT panmat-targets
    FILE panmat-targets.cmake
    NAMESPACE panmat::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/panmat
)

include(CMakePackageConfigHelpers)

include(GNUInstallDirs)

set(INCLUDE_INSTALL_DIR 
    CACHE PATH "Location of header files" )

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/panmat-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/panmat-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/panmat
	PATH_VARS INCLUDE_INSTALL_DIR
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/panmat-config.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/panmat"
)