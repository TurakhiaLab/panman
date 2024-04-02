cmake_minimum_required (VERSION 3.27) 

project(panmat
	VERSION 0.1.0
	DESCRIPTION "Pangenome mutation-annotated trees"
	LANGUAGES CXX)

#  To consume the panmat library, use:
#
# (a) find_package(panmat [HINTS path/to/panmat])
#     after installing this package with `make install`, or
#
# (b) FetchContent_Declare(panmat
#	    GIT_REPOSITORY https://github.com/TurakhiaLab/pagenome-mat)
# then: FetchContent_MakeAvailable(panmat)

set(CMAKE_CXX_STANDARD 17) # could bump to 20?

set(BUILD_TESTING OFF)

include(FetchContent)

## Main targets

# pan-mat library
add_library(panmat src/PangenomeMAT.cpp src/PangenomeMAT.hpp)

# application binary 
add_executable(panmat-utils src/main.cpp)

target_include_directories(panmat-utils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

find_package(ZLIB REQUIRED GLOBAL)
message(" ✔︎ found local zlib")

# jsoncpp
find_package(jsoncpp QUIET GLOBAL)
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
		OVERRIDE_FIND_PACKAGE TRUE
	)
	FetchContent_MakeAvailable(jsoncpp)
	message(" ✔︎ using fetched jsoncpp")
	list(APPEND REMOTE_DEPS jsoncpp_static)
else()
	list(APPEND LOCAL_DEPS jsoncpp_static)
	message(" ✔︎ found local jsoncpp")
endif()

# protobuf
find_package(Protobuf 5.24.3 QUIET CONFIG)
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
		OVERRIDE_FIND_PACKAGE TRUE
	)
	FetchContent_MakeAvailable(protobuf)
	include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
endif()
protobuf_generate(TARGET panmat
                  PROTOS mutation_annotation_test_proto3_optional_new.proto
                  LANGUAGE cpp
                  IMPORT_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if (Protobuf_USE_LOCAL)
	list(APPEND LOCAL_DEPS ${Protobuf_LIBRARIES} ${PROTOBUF_LIBRARY})
	message(" ✔︎ found local protobuf")
else()
	list(APPEND REMOTE_DEPS ${Protobuf_LIBRARIES} ${PROTOBUF_LIBRARY})
	message(" ✔︎ using fetched protobuf")
endif()
message("Protobuf_LIBRARIES: ${Protobuf_LIBRARIES}")
# boost
find_package(Boost 1.75.0 QUIET COMPONENTS system filesystem program_options iostreams date_time GLOBAL)
if (NOT ${Boost_FOUND})
	set(BOOST_INCLUDE_LIBRARIES system filesystem program_options iostreams date_time)
	set(BOOST_ENABLE_CMAKE ON)
	FetchContent_Declare(
		Boost
		GIT_REPOSITORY https://github.com/boostorg/boost
		GIT_TAG ad09f667e61e18f5c31590941e748ac38e5a81bf # 1.84.0
		GIT_SHALLOW TRUE
		OVERRIDE_FIND_PACKAGE TRUE
	)
	list(APPEND REMOTE_DEPS Boost::program_options Boost::filesystem Boost::iostreams Boost::date_time Boost::system)
	FetchContent_MakeAvailable(Boost)
	message(" ✔︎ using fetched boost")
else()
	list(APPEND LOCAL_DEPS Boost::program_options Boost::filesystem Boost::iostreams Boost::date_time Boost::system)
	message(" ✔︎ found local boost")
endif()

# spoa
find_package(spoa QUIET GLOBAL)
if (NOT ${spoa_FOUND})
	FetchContent_Declare(
		spoa
		GIT_REPOSITORY https://github.com/rvaser/spoa
		GIT_TAG f39e67f229168cc305e2d5cc8765cd30c8e308d5 # v4.1.4
		GIT_SHALLOW TRUE
		OVERRIDE_FIND_PACKAGE TRUE
	)
	list(APPEND REMOTE_DEPS spoa)
	FetchContent_MakeAvailable(spoa)
	message(" ✔︎ using fetched spoa")
else()
	list(APPEND LOCAL_DEPS spoa)
	message(" ✔︎ found local spoa")
endif()

# tbb
find_package(TBB 2021.11.0 QUIET GLOBAL)
if(NOT ${TBB_FOUND})
	set(TBB_TEST OFF CACHE INTERNAL "")
	set(TBB_FIND_PACKAGE OFF CACHE INTERNAL "")
	set(BUILD_SHARED_LIBS ON CACHE INTERNAL "")
	set(TBB_EXAMPLES OFF CACHE INTERNAL "")
	set(TBB_INCLUDE_DIRS ${FETCHCONTENT_BASE_DIR}/tbb-src/include CACHE PATH "" FORCE)
	FetchContent_Declare(
		TBB
		GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
		GIT_TAG 8b829acc65569019edb896c5150d427f288e8aba
		GIT_SHALLOW TRUE
	)
	include_directories(${TBB_INCLUDE_DIRS})
	list(APPEND REMOTE_DEPS TBB::tbb)
	FetchContent_MakeAvailable(TBB)
	message(" ✔︎ using fetched tbb")
else()
	list(APPEND LOCAL_DEPS TBB::tbb)
	message(" ✔︎ found local tbb")
endif()


# public header
set_target_properties(panmat PROPERTIES PUBLIC_HEADER src/PangenomeMAT.hpp)

## Dependencies

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE BOOL "" FORCE)

## Linking

target_link_libraries(panmat PUBLIC
	Boost::program_options Boost::filesystem Boost::iostreams Boost::date_time Boost::system
	protobuf::libprotobuf
	TBB::tbb
	spoa
	jsoncpp_static
	ZLIB::ZLIB
)

target_link_libraries(panmat-utils PUBLIC panmat)

## Exporting
#   - `make` will create build/panmat-config.cmake (if make is run)
#   - {install path}/cmake/panmat-config.cmake (if make install is run)

install(TARGETS panmat-utils
	EXPORT panmat-targets
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include
)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/panmat-config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/panmat-config.cmake"
    INSTALL_DESTINATION lib/cmake/panmat
)

install(EXPORT panmat-targets
    FILE panmat-targets.cmake
    NAMESPACE panmat::
    DESTINATION lib/cmake/panmat
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/panmat-config.cmake"
    DESTINATION lib/cmake/panmat
)