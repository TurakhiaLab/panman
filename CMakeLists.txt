cmake_minimum_required (VERSION 3.27) 

project(panmat)

# Installs pan-MAT library to "include/PangenomeMAT.hpp"
# and panmat-utils executable to "bin/panmat-utils"

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# main target
add_executable(panmat-utils
	src/main.cpp
	src/PangenomeMAT.cpp
	src/PangenomeMAT.hpp)

# zlib
find_package(ZLIB REQUIRED)
message(" ✔︎ found local zlib")


# jsoncpp
find_package(jsoncpp QUIET)
if (NOT ${jsoncpp_FOUND})
	FetchContent_Declare(
		jsoncpp
		GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp
		OVERRIDE_FIND_PACKAGE TRUE
	)
	FetchContent_MakeAvailable(jsoncpp)
	find_package(jsoncpp REQUIRED) # just prints version due to override
	message(" ✔︎ using fetched jsoncpp")
else()
	message(" ✔︎ found local jsoncpp")
endif()

# protobuf
find_package(Protobuf 5.24.3 QUIET NAMES Protobuf protobuf)
if (Protobuf_FOUND)
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
if(NOT Protobuf_USE_LOCAL)
	set(ABSL_PROPAGATE_CXX_STD ON)
	set(protobuf_BUILD_TESTS OFF)
	set(protobuf_BUILD_CONFORMANCE OFF)
	set(protobuf_BUILD_EXAMPLES OFF)
	set(protobuf_BUILD_PROTOBUF_BINARIES ON)
	set(protobuf_BUILD_SHARED_LIBS OFF)
	set(protobuf_MSVC_STATIC_RUNTIME OFF)
    FetchContent_Declare(
        protobuf
        GIT_REPOSITORY https://github.com/protocolbuffers/protobuf
        GIT_TAG d6511091a0cab1ad13f676a02676ad2a0e5eb9ae # Specify the version you need 
		GIT_SHALLOW TRUE
		OVERRIDE_FIND_PACKAGE TRUE
	)
	list(APPEND CMAKE_PREFIX_PATH "${protobuf_BINARY_DIR}")
	FetchContent_MakeAvailable(protobuf)
	find_package(Protobuf REQUIRED) # just prints version due to override
	include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
endif()
protobuf_generate(TARGET panmat-utils
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
find_package(Boost QUIET COMPONENTS system filesystem program_options iostreams date_time)

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
	FetchContent_MakeAvailable(Boost)
	find_package(Boost REQUIRED) # just prints version due to override
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
		OVERRIDE_FIND_PACKAGE TRUE
	)
	FetchContent_MakeAvailable(spoa)
	find_package(spoa REQUIRED) # just prints version due to override
	message(" ✔︎ using fetched spoa")
else()
	message(" ✔︎ found local spoa")
endif()

# tbb
find_package(TBB 2021.10.0 QUIET)
if(NOT ${TBB_FOUND})
	set(FETCHCONTENT_QUIET FALSE)
	set(_GLIBCXX_USE_CXX11_ABI FALSE)
	FetchContent_Declare(
		TBB
		GIT_REPOSITORY https://github.com/oneapi-src/oneTBB/
		GIT_TAG 8b829acc65569019edb896c5150d427f288e8aba # v2021.11.0
		GIT_SHALLOW TRUE
		OVERRIDE_FIND_PACKAGE TRUE
	)
	FetchContent_MakeAvailable(TBB)
	find_package(TBB REQUIRED) # just prints version due to override
	message(" ✔︎ using fetched tbb")
else()
	message(" ✔︎ found local tbb")
endif()

# targets
set(DEP_TARGETS
	Boost::system Boost::filesystem Boost::program_options Boost::iostreams Boost::date_time
	TBB::tbb
	protobuf::libprotobuf
	spoa
	jsoncpp_static
	z)

target_link_libraries(panmat-utils PUBLIC ${DEP_TARGETS})

set_target_properties(panmat-utils PROPERTIES PUBLIC_HEADER src/PangenomeMAT.hpp)

export(TARGETS panmat-utils 
  FILE "${CMAKE_CURRENT_BINARY_DIR}/panmat-utils-targets.cmake"
)	
install(TARGETS panmat-utils
	EXPORT panmat-utils-targets
 	ARCHIVE DESTINATION lib
 	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
 	PUBLIC_HEADER DESTINATION include
 )