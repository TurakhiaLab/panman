include(FetchContent)
include(ExternalProject)

set(BUILD_TESTING OFF)
set(BUILD_STATIC_LIBS ON)
set(POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS OFF)

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
FetchContent_Declare(
	jsoncpp
	GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp
	FIND_PACKAGE_ARGS)

# protobuf
set(protobuf_BUILD_TESTS OFF CACHE INTERNAL "")
set(protobuf_BUILD_EXPORT ON CACHE INTERNAL "")
set(protobuf_MSVC_STATIC_RUNTIME OFF CACHE INTERNAL "")
set(protobuf_BUILD_PROTOBUF_BINARIES ON CACHE INTERNAL "")
set(protobuf_BUILD_LIBPROTOC OFF CACHE INTERNAL "")
set(protobuf_BUILD_SHARED_LIBS OFF CACHE INTERNAL "")
set(ABSL_BUILD_DLL OFF CACHE BOOL "" FORCE)
set(ABSL_PROPAGATE_CXX_STD ON CACHE INTERNAL "")
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


set(BUILD_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)

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

include(${protobuf_SOURCE_DIR}/cmake/protobuf-generate.cmake)
set(protobuf_INCLUDE_DIRS ${protobuf_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)

FetchContent_MakeAvailable(Boost TBB ZLIB spoa protobuf jsoncpp)

