cmake_minimum_required (VERSION 3.8) 

project(panmanUtils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

if(DEFINED Protobuf_PATH) 
	find_package(Protobuf REQUIRED HINTS ${Protobuf_PATH}) 
else()
	find_package(Protobuf REQUIRED) 
endif()

# Print version messages
if(Protobuf_FOUND)
	message(STATUS "Using Protocol Buffers ${Protobuf_VERSION}")
endif()
include_directories(${Protobuf_INCLUDE_DIRS})

# Include TBB
# include(${TBB_DIR}/cmake/TBBBuild.cmake)
# if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm")
# 	set(ARCH "arm64")
# elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "armv4t")
# 	set(ARCH "arm64")
# elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
# 	set(ARCH "x86")
# elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
# 	set(ARCH "x86_64")
# endif()

# tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1 arch=${ARCH})
# find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview HINTS ${TBB_DIR})

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(BUILD_STATIC_LIBS ON)
set(BUILD_SHARED_LIBS ON)
set(JSONCPP_WITH_TESTS OFF)
set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF)
set(JSONCPP_WITH_CMAKE_PACKAGE ON)
set(ARCHIVE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
find_package(PkgConfig REQUIRED)

FetchContent_Declare(
	jsoncpp
	GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
	GIT_TAG	1.9.5
	GIT_SHALLOW	TRUE
)

FetchContent_MakeAvailable(jsoncpp)

# pkg_check_modules(JSONCPP jsoncpp_static IMPORTED_TARGET GLOBAL)

file(GLOB PANMAT_SRCS "src/panmanUtils.cpp" "src/panman.cpp" "src/panmanUtils.hpp" "src/panman.hpp")


	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		panman.proto)
	
	add_executable(panmanUtils
		${PANMAT_SRCS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)

target_include_directories(panmanUtils PUBLIC ${TBB_DIR}/include)

TARGET_COMPILE_OPTIONS(panmanUtils PUBLIC -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

include(FetchContent)
FetchContent_Declare(
	TBB
	GIT_REPOSITORY https://github.com/wjakob/tbb
	GIT_TAG ddbe45cd3ad89df9a84cd77013d5898fc48b8e89
	GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(TBB)
TARGET_LINK_LIBRARIES(panmanUtils PRIVATE  ${Boost_LIBRARIES} ${TBB_LIBRARIES} ${Protobuf_LIBRARIES} jsoncpp_static stdc++) #${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
target_include_directories(panmanUtils PUBLIC "${PROJECT_BINARY_DIR}")


add_library(panman SHARED ${PANMAT_SRCS} ${PROTO_SRCS} ${PROTO_HDRS})
target_include_directories(panman PUBLIC "${PROJECT_BINARY_DIR}")
target_include_directories(panman PUBLIC ${TBB_DIR}/include)

target_link_libraries(panman PRIVATE ${Boost_LIBRARIES} ${TBB_LIBRARIES} ${Protobuf_LIBRARIES} jsoncpp_static stdc++)
