cmake_minimum_required (VERSION 3.8) 

project(panmat-utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
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

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include(${TBB_DIR}/cmake/TBBBuild.cmake)
include("/home/AD.UCSD.EDU/hmotwani/vcpkg/scripts/buildsystems/vcpkg.cmake")
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)
find_package(jsoncpp CONFIG REQUIRED)
# target_link_libraries(main PRIVATE JsonCpp::JsonCpp)

#find_package(ZLIB REQUIRED)


file(GLOB PANMAT_SRCS "src/main.cpp" "src/PangenomeMAT.cpp" "src/PangenomeMAT.hpp")

if(DEFINED Protobuf_PATH)
	add_executable(panmat-utils
		${PANMAT_SRCS}
		)

	protobuf_generate(
		LANGUAGE cpp
		TARGET panmat-utils 
		PROTOS mutation_annotation_test_proto3_optional.proto mutation_annotation_test_proto3_optional_new.proto vg.proto)

	# protobuf_generate(
	# 	LANGUAGE cpp
	# 	TARGET panmat-utils 
	# 	PROTOS mutation_annotation_test_proto3_optional_new.proto)


else()
	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		mutation_annotation_test_proto3_optional.proto mutation_annotation_test_proto3_optional_new.proto vg.proto)
	
	# protobuf_generate_cpp(
	# 	PROTO_SRCS PROTO_HDRS
	# 	mutation_annotation_test_proto3_optional_new.proto)

	add_executable(panmat-utils
		${PANMAT_SRCS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)

endif()

TARGET_COMPILE_OPTIONS(panmat-utils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)

add_subdirectory(spoa EXCLUDE_FROM_ALL)

# find_package(spoa)
TARGET_LINK_LIBRARIES(panmat-utils PRIVATE spoa::spoa stdc++  JsonCpp::JsonCpp ${Boost_LIBRARIES} -lboost_iostreams ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES}) #${Boost_LIBRARIES}  ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
target_include_directories(panmat-utils PUBLIC "${PROJECT_BINARY_DIR}")
