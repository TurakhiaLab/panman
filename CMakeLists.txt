cmake_minimum_required (VERSION 3.27) 

project(panmat
	VERSION 0.1.0
	DESCRIPTION "Pangenome mutation-annotated trees"
	LANGUAGES CXX)

# Installs panmat library to "include/PangenomeMAT.hpp"
# and panmat executable to "bin/panmat-utils"

# To consume the panmat library, use 
# (a) find_package(panmat [HINTS path/to/panmat]) or
# (b) FetchContent_Declare(panmat GIT_REPOSITORY https://github.com/TurakhiaLab/pagenome-mat)
# 	  then FetchContent_MakeAvailable(panmat)

set(CMAKE_CXX_STANDARD 17) # could bump to 20

set(BUILD_TESTING OFF)

include(FetchContent)

## Main targets

add_library(panmat)

add_executable(panmat-utils
	src/main.cpp
	src/PangenomeMAT.cpp
	src/PangenomeMAT.hpp)

## Dependencies

# zlib
find_package(ZLIB REQUIRED)
message(" ✔︎ found local zlib")

# jsoncpp
find_package(jsoncpp QUIET)
if (NOT ${jsoncpp_FOUND})
	set(JSONCPP_WITH_TESTS OFF CACHE INTERNAL "")
	set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE INTERNAL "")
	set(JSONCPP_WITH_PKGCONFIG_SUPPORTOFF CACHE INTERNAL "")
	set(JSONCPP_WITH_STRICT_ISO CACHE INTERNAL "")
	set(BUILD_STATIC_LIBS ON CACHE INTERNAL "")
	set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "")
	set(BUILD_OBJECT_LIBS OFF CACHE INTERNAL "")
	FetchContent_Declare(
		jsoncpp
		GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp
	)
	FetchContent_MakeAvailable(jsoncpp)
	message(" ✔︎ using fetched jsoncpp")
else()
	message(" ✔︎ found local jsoncpp")
endif()

# protobuf
find_package(Protobuf 5.24.3 QUIET NAMES Protobuf protobuf)
if (Protobuf_FOUND AND Protobuf_PROTOC_EXECUTABLE)
	execute_process(
		COMMAND ${Protobuf_PROTOC_EXECUTABLE} --version
		OUTPUT_VARIABLE PROTOC_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	# protoc compiler version should match library
	string(SUBSTRING ${PROTOC_VERSION} 10 -1 PROTOC_VERSION)
	string(SUBSTRING ${Protobuf_VERSION} 2 -1 Protobuf_VERSION_MINOR)
	if (Protobuf_VERSION_MINOR EQUAL PROTOC_VERSION)
		set(Protobuf_USE_LOCAL ON)
	endif()
endif()
if(NOT Protobuf_FOUND AND NOT Protobuf_USE_LOCAL)
	set(ABSL_PROPAGATE_CXX_STD ON)
	set(protobuf_BUILD_TESTS OFF CACHE INTERNAL "")
	set(protobuf_BUILD_PROTOBUF_BINARIES ON CACHE INTERNAL "")
	set(protobuf_BUILD_PROTOC_BINARIES ON CACHE INTERNAL "")
	set(protobuf_BUILD_LIBPROTOC OFF CACHE INTERNAL "")
	set(protobuf_BUILD_SHARED_LIBS OFF CACHE INTERNAL "")
	FetchContent_Declare(
        protobuf
        GIT_REPOSITORY https://github.com/protocolbuffers/protobuf
        GIT_TAG d6511091a0cab1ad13f676a02676ad2a0e5eb9ae # Specify the version you need 
		GIT_SHALLOW TRUE
	)
	list(APPEND CMAKE_PREFIX_PATH "${protobuf_BINARY_DIR}")
	FetchContent_MakeAvailable(protobuf)
	include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
endif()
protobuf_generate(TARGET panmat
                  PROTOS mutation_annotation_test_proto3_optional_new.proto
                  LANGUAGE cpp
                  IMPORT_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
if (Protobuf_USE_LOCAL)
	message(" ✔︎ found local protobuf")
else()
	message(" ✔︎ using fetched protobuf")
endif()

# boost
find_package(Boost 1.75.0 QUIET COMPONENTS system filesystem program_options iostreams date_time)
if (NOT ${Boost_FOUND})
	set(BOOST_INCLUDE_LIBRARIES system filesystem program_options iostreams date_time)
	set(BOOST_ENABLE_CMAKE ON)
	FetchContent_Declare(
		Boost
		GIT_REPOSITORY https://github.com/boostorg/boost
		GIT_TAG ad09f667e61e18f5c31590941e748ac38e5a81bf # 1.84.0
		GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(Boost)
	message(" ✔︎ using fetched boost")
else()
	message(" ✔︎ found local boost")
endif()

# spoa
find_package(spoa QUIET)
if (NOT ${spoa_FOUND})
	FetchContent_Declare(
		spoa
		GIT_REPOSITORY https://github.com/rvaser/spoa
		GIT_TAG f39e67f229168cc305e2d5cc8765cd30c8e308d5 # v4.1.4
		GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(spoa)
	message(" ✔︎ using fetched spoa")
else()
	message(" ✔︎ found local spoa")
endif()

# tbb
find_package(TBB 2021.12.0 QUIET)
if(NOT ${TBB_FOUND})
	set(TBB_TEST OFF CACHE INTERNAL "")
	set(TBB_BENCH OFF CACHE INTERNAL "")
	set(TBB_BUILD ON CACHE INTERNAL "")
	set(TBB_FIND_PACKAGE OFF CACHE INTERNAL "")
	set(BUILD_SHARED_LIBS ON CACHE INTERNAL "")
	set(TBB_EXAMPLES OFF CACHE INTERNAL "")
	set(TBB_INCLUDE_DIRS ${FETCHCONTENT_BASE_DIR}/tbb-src/include CACHE PATH "" FORCE)
	block()
		set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE BOOL "" FORCE)
		FetchContent_Declare(
			TBB
			GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
			GIT_TAG 8b829acc65569019edb896c5150d427f288e8aba
			GIT_SHALLOW TRUE
		)
		FetchContent_MakeAvailable(TBB)
	endblock()
	include_directories(${TBB_INCLUDE_DIRS})
	message(" ✔︎ using fetched tbb")
else()
	message(" ✔︎ found local tbb")
endif()

## Linking

set(DEP_TARGETS
	Boost::system Boost::filesystem Boost::program_options Boost::iostreams Boost::date_time
	# tbb cmake integration is weird, this should cover all names it could have
	$<TARGET_NAME_IF_EXISTS:TBB> $<TARGET_NAME_IF_EXISTS:tbb> $<TARGET_NAME_IF_EXISTS:TBB::tbb>
	$<TARGET_NAME_IF_EXISTS:tbbmalloc> $<TARGET_NAME_IF_EXISTS:TBB::tbbmalloc>
	$<TARGET_NAME_IF_EXISTS:tbbmalloc_proxy> $<TARGET_NAME_IF_EXISTS:TBB::tbbmalloc_proxy>
	protobuf::libprotobuf
	spoa
	jsoncpp_static
	z)

target_link_libraries(panmat PUBLIC ${DEP_TARGETS})

target_link_libraries(panmat-utils PUBLIC panmat)

## Exporting

# Set up config files
# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(
# 	${CMAKE_CURRENT_BINARY_DIR}/panmat-config-version.cmake
# 	VERSION ${PROJECT_VERSION}
# 	COMPATIBILITY AnyNewerVersion
# )

# # Include directories
# target_include_directories(panmat PUBLIC
# 	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
# 	$<INSTALL_INTERFACE:include>
# )

# # Install targets (for make install)
# install(TARGETS panmat spoa
# 	EXPORT panmat-targets
# 	LIBRARY DESTINATION lib
# 	ARCHIVE DESTINATION lib
# 	INCLUDES DESTINATION include
# )
# install(TARGETS panmat-utils
# 	EXPORT panmat-targets
# 	RUNTIME DESTINATION bin
# )

# install(EXPORT panmat-targets
# 	FILE panmat-targets.cmake
# 	NAMESPACE panmat::
# 	DESTINATION lib/cmake/panmat
# )

# configure_file(panmat-config.cmake.in panmat-config.cmake @ONLY)

# install(FILES "${CMAKE_CURRENT_BINARY_DIR}/panmat-config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/panmat-config-version.cmake"
# 	DESTINATION lib/cmake/panmat
# )