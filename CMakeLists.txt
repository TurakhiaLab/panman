cmake_minimum_required (VERSION 3.20)

cmake_policy(SET CMP0135 NEW)

project(panmat
	VERSION 0.1.0
	DESCRIPTION "Pangenome mutation-annotated trees"
	LANGUAGES CXX)

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
#

set(CMAKE_CXX_STANDARD 17) # could bump to 20?

set(BUILD_TESTING OFF)

include(FetchContent)

## Main targets ##

# pan-mat library
add_library(panmat src/PangenomeMAT.cpp src/PangenomeMAT.hpp)
# application binary 
add_executable(panmat-utils src/main.cpp src/PangenomeMAT.cpp src/PangenomeMAT.hpp)
## Dependencies ##

# zlib
FetchContent_Declare(ZLIB
    URL https://www.zlib.net/zlib-1.2.13.tar.gz
        https://www.zlib.net/fossils/zlib-1.2.13.tar.gz
        https://github.com/madler/zlib/releases/download/v1.2.13/zlib-1.2.13.tar.gz
    URL_HASH MD5=9b8aa094c4e5765dabf4da391f00d15c
	FIND_PACKAGE_ARGS
	DOWNLOAD_EXTRACT_TIMESTAMP
)

# jsoncpp
set(JSONCPP_WITH_TESTS OFF CACHE INTERNAL "")
set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE INTERNAL "")
set(JSONCPP_WITH_PKGCONFIG_SUPPORTOFF CACHE INTERNAL "")
set(JSONCPP_WITH_STRICT_ISO CACHE INTERNAL "")
# set(BUILD_STATIC_LIBS OFF CACHE INTERNAL "")
# set(BUILD_SHARED_LIBS ON CACHE INTERNAL "")
# set(BUILD_OBJECT_LIBS OFF CACHE INTERNAL "")
FetchContent_Declare(
	jsoncpp
	GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp
	FIND_PACKAGE_ARGS)

# protobuf
FetchContent_Declare(
	protobuf
	GIT_REPOSITORY https://github.com/protocolbuffers/protobuf
	GIT_TAG 2434ef2adf0c74149b9d547ac5fb545a1ff8b6b5 # Specify the version you need 
	GIT_SHALLOW TRUE
	FIND_PACKAGE_ARGS 5.25.3 CONFIG
)

# boost
set(BOOST_INCLUDE_LIBRARIES system filesystem program_options iostreams date_time)
set(BOOST_ENABLE_CMAKE ON)
FetchContent_Declare(
	Boost
	GIT_REPOSITORY https://github.com/boostorg/boost
	GIT_TAG ad09f667e61e18f5c31590941e748ac38e5a81bf # 1.84.0
	GIT_SHALLOW TRUE
	FIND_PACKAGE_ARGS 1.75.0 COMPONENTS system filesystem program_options iostreams date_time
)

# spoa
FetchContent_Declare(
	spoa
	GIT_REPOSITORY https://github.com/rvaser/spoa
	GIT_TAG f39e67f229168cc305e2d5cc8765cd30c8e308d5 # v4.1.4
	GIT_SHALLOW TRUE
	FIND_PACKAGE_ARGS
)

# tbb
set(TBB_TEST OFF CACHE INTERNAL "")	
set(TBB_FIND_PACKAGE OFF CACHE INTERNAL "")
set(TBB_EXAMPLES OFF CACHE INTERNAL "")
set(TBB_INCLUDE_DIRS ${FETCHCONTENT_BASE_DIR}/tbb-src/include CACHE PATH "" FORCE)
FetchContent_Declare(
	TBB
	GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
	GIT_TAG 8b829acc65569019edb896c5150d427f288e8aba
	GIT_SHALLOW TRUE
	FIND_PACKAGE_ARGS 2021.11.0
)
include_directories(${TBB_INCLUDE_DIRS})

target_include_directories(panmat PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Run protoc compiler
set(panmat_SRCS mutation_annotation_test_proto3_optional_new.pb.cc)
set(panmat_PROTOS mutation_annotation_test_proto3_optional_new.proto)

target_sources(panmat PRIVATE ${panmat_SRCS} ${panmat_PROTOS})

#set(CMAKE_INCLUDE_CURRENT_DIR TRUE CACHE INTERNAL "")
set(protobuf_BUILD_TESTS OFF CACHE INTERNAL "")
set(protobuf_BUILD_EXPORT ON CACHE INTERNAL "")
set(protobuf_MSVC_STATIC_RUNTIME OFF CACHE INTERNAL "")
set(protobuf_BUILD_PROTOBUF_BINARIES ON CACHE INTERNAL "")
set(protobuf_BUILD_LIBPROTOC OFF CACHE INTERNAL "")
set(protobuf_BUILD_SHARED_LIBS OFF CACHE INTERNAL "")
set(POSITION_INDEPENDENT_CODE ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(ABSL_BUILD_DLL OFF CACHE BOOL "" FORCE)
set(ABSL_PROPAGATE_CXX_STD ON CACHE INTERNAL "")

FetchContent_MakeAvailable(jsoncpp protobuf)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
set(protobuf_INCLUDE_DIRS ${protobuf_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)

include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
#list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)

protobuf_generate(TARGET panmat)

FetchContent_MakeAvailable(Boost TBB ZLIB spoa)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE BOOL "" FORCE)

## Linking ##

target_link_libraries(panmat PUBLIC
	Boost::program_options Boost::filesystem Boost::iostreams Boost::date_time Boost::system
	protobuf::libprotobuf
	TBB::tbb
	$<TARGET_NAME_IF_EXISTS:jsoncpp_static> $<TARGET_NAME_IF_EXISTS:jsoncpp_lib>
	spoa
	ZLIB::ZLIB
)

target_link_libraries(panmat-utils PRIVATE panmat)


## Exporting
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
