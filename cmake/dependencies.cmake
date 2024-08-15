include(FetchContent)
include(ExternalProject)

set(FETCHCONTENT_QUIET OFF)
set(BUILD_TESTING OFF)
set(BUILD_STATIC_LIBS ON)
set(BUILD_SHARED_LIBS OFF)

set(CAN_FIND_PACKAGE OFF)
if ("${CMAKE_VERSION}" VERSION_LESS "3.24.0")
    set(CAN_FIND_PACKAGE ON)
endif()

#write generator expressions to file

# jsoncpp
set(JSONCPP_WITH_TESTS OFF CACHE INTERNAL "")
set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE INTERNAL "")
set(JSONCPP_WITH_PKGCONFIG_SUPPORTOFF CACHE INTERNAL "")
set(JSONCPP_WITH_STRICT_ISO CACHE INTERNAL "")
FetchContent_Declare(
	jsoncpp
	GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp
    GIT_SHALLOW TRUE
)

# protobuf
set(protobuf_BUILD_TESTS OFF CACHE INTERNAL "")
set(protobuf_BUILD_EXPORT OFF CACHE INTERNAL "")
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
)

# spoa
FetchContent_Declare(
	spoa
	GIT_REPOSITORY https://github.com/rvaser/spoa
	GIT_TAG f39e67f229168cc305e2d5cc8765cd30c8e308d5 # v4.1.4
	GIT_SHALLOW TRUE
)

# shared libraries
set(BUILD_STATIC_LIBS OFF)
set(POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS ON)

# boost
set(BOOST_INCLUDE_LIBRARIES system filesystem iostreams program_options)
set(BOOST_ENABLE_CMAKE ON)
FetchContent_Declare(
	Boost
	URL https://github.com/boostorg/boost/releases/download/boost-1.85.0.beta1/boost-1.85.0.beta1-cmake.tar.gz
    URL_HASH SHA256=80014af00b69f6d6bb20eb1ae659cf24e6c439537e88e99ccc03fc4bae7ffe84
    GIT_SHALLOW TRUE
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
)
include_directories(${TBB_INCLUDE_DIRS})

set(protobuf_INCLUDE_DIRS ${protobuf_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
