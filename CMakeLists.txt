cmake_minimum_required (VERSION 3.27) 

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
project(panmat-utils)

include(FetchContent)


# zlib
find_package(ZLIB REQUIRED)
message(" ✔️ zlib ${ZLIB_VERSION_STRING}")

# protobuf
set(ABSL_PROPAGATE_CXX_STD ON)
find_package(Protobuf QUIET)
if ("${Protobuf_FOUND}")
	execute_process(
		COMMAND ${Protobuf_PROTOC_EXECUTABLE} --version
		OUTPUT_VARIABLE PROTOC_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	string(SUBSTRING ${PROTOC_VERSION} 10 -1 PROTOC_VERSION)
	string(SUBSTRING ${Protobuf_VERSION} 2 -1 Protobuf_VERSION_MINOR)
	if (Protobuf_VERSION_MINOR EQUAL PROTOC_VERSION)
		set(Protobuf_USE_LOCAL ON)
	endif()	
endif()
if(NOT ${Protobuf_FOUND} OR NOT ${Protobuf_USE_LOCAL})
	message(" ↯ (using fetched protobuf)")
	FetchContent_Declare(
		Protobuf
		GIT_REPOSITORY "https://github.com/protocolbuffers/protobuf.git"
		GIT_TAG d6511091a0cab1ad13f676a02676ad2a0e5eb9ae # v5.26.0
		GIT_SHALLOW TRUE
	  )
	  set(TMP ${CMAKE_EXE_LINKER_FLAGS})
	  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -w")
		  FetchContent_MakeAvailable(Protobuf)
	  set(CMAKE_EXE_LINKER_FLAGS ${TMP})
	  #FetchContent_GetProperties(Protobuf SOURCE_DIR Protobuf_SOURCE_DIR)
	  FetchContent_GetProperties(Protobuf INCLUDE_DIR Protobuf_INCLUDE_DIRS)
	  set(Protobuf_VERSION 26.0)
endif()
include_directories(${Protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
message(" ✔️ protobuf ${Protobuf_VERSION}")

protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS mutation_annotation_test_proto3_optional_new.proto)

# boost
find_package(Boost QUIET COMPONENTS system filesystem program_options iostreams date_time)

if (NOT ${Boost_FOUND})
	message(" ↯ (using fetched boost)")
	set(BOOST_INCLUDE_LIBRARIES system filesystem program_options iostreams date_time)
	set(BOOST_ENABLE_CMAKE ON)
	FetchContent_Declare(
		Boost
		GIT_REPOSITORY https://github.com/boostorg/boost.git
		GIT_TAG ad09f667e61e18f5c31590941e748ac38e5a81bf # 1.84.0
		GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(Boost)
	set(Boost_VERSION 1.84.0)
endif()
include_directories(${Boost_INCLUDE_DIRS})
message(" ✔️ boost ${Boost_VERSION}")


# tbb
find_package(TBB 2021)

if(NOT ${TBB_FOUND})
	message(" ↯ (using fetched tbb)")
	set(TBB_BUILD_SHARED OFF)
	set(TBB_BUILD_STATIC ON)
	set(TBB_BUILD_TESTS OFF)
	set(TBB_TEST OFF)
	set(TBB_BUILD_TBBMALLOC ON)
	FetchContent_Declare(
		TBB
		GIT_REPOSITORY https://github.com/oneapi-src/oneTBB/
		GIT_TAG 8b829acc65569019edb896c5150d427f288e8aba # v2021.11.0
		GIT_SHALLOW TRUE
	)
	cmake_policy(SET CMP0063 NEW) # suppress tbb warning
	FetchContent_MakeAvailable(TBB)
	set(TBB_VERSION 2021.11.0)
	FetchContent_GetProperties(TBB INCLUDE_DIR TBB_INCLUDE_DIRS)
endif()
include_directories(${TBB_INCLUDE_DIRS})

message(" ✔️ tbb ${TBB_VERSION}")

add_executable(panmat-utils
	${PROTO_SRCS}
	${PROTO_HDRS}
	src/main.cpp
	src/PangenomeMAT.cpp
	src/PangenomeMAT.hpp
)
add_library(panmat
	${PROTO_SRCS}
	${PROTO_HDRS}
	src/PangenomeMAT.cpp
	src/PangenomeMAT.hpp
)
# jsoncpp
find_package(jsoncpp CONFIG QUIET)
if (NOT ${jsoncpp_FOUND})
	message(" ↯ (using fetched jsoncpp)")
	FetchContent_Declare(
		jsoncpp
		GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp
		GIT_TAG 5defb4ed1a4293b8e2bf641e16b156fb9de498cc # 1.9.5
	)
	FetchContent_MakeAvailable(jsoncpp)
	set(jsoncpp_VERSION 1.9.5)
endif()
include_directories(${jsoncpp_INCLUDE_DIRS})
message(" ✔️ jsoncpp ${jsoncpp_VERSION}")

# spoa
find_package(spoa QUIET)
if (NOT ${spoa_FOUND})
	message(" ↯ (using fetched spoa)")
	FetchContent_Declare(
		spoa
		GIT_REPOSITORY https://github.com/rvaser/spoa
		GIT_TAG f39e67f229168cc305e2d5cc8765cd30c8e308d5) # v4.1.4
	FetchContent_MakeAvailable(spoa)
	set(spoa_VERSION 4.1.4)
endif()
include_directories(${spoa_INCLUDE_DIRS})
message(" ✔️ spoa ${spoa_VERSION}")

target_link_libraries(panmat-utils PUBLIC
	protobuf::libprotobuf
	spoa
	jsoncpp
	Boost::system Boost::filesystem Boost::program_options Boost::iostreams Boost::date_time
	TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy
	z
)
set_target_properties(panmat PROPERTIES PUBLIC_HEADER src/PangenomeMAT.hpp)
export(TARGETS panmat FILE panmatConfig.cmake)
install(TARGETS panmat
 	ARCHIVE DESTINATION "lib"
 	LIBRARY DESTINATION "lib"
 	PUBLIC_HEADER DESTINATION "include"
 )
