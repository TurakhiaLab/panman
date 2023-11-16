cmake_minimum_required (VERSION 3.16) 

project(panmat-utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++17 -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_INCLUDE_CURRENT_DIR ON)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")



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

tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)

find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

#find_package(ZLIB)



add_custom_target(libminimap2.a
   COMMAND make libminimap2.a
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/minimap2_src/
)

add_custom_target(clean-minimap
    COMMAND make clean
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/minimap2_src/
)



file(GLOB PANMAT_SRCS "src/align_from_seeds.c" "src/minimap2_src/libminimap2.a" "src/*.cpp" "src/*.hpp")



if(DEFINED Protobuf_PATH)
	add_executable(panmat-utils
		${PANMAT_SRCS}
		)

	protobuf_generate(
		LANGUAGE cpp
		TARGET panmat-utils 
		PROTOS mutation_annotation_test_proto3_optional.proto mutation_annotation_test_proto3_optional_new.proto)

	# protobuf_generate(
	# 	LANGUAGE cpp
	# 	TARGET panmat-utils 
	# 	PROTOS mutation_annotation_test_proto3_optional_new.proto)


else()
	protobuf_generate_cpp(
		PROTO_SRCS PROTO_HDRS
		mutation_annotation_test_proto3_optional.proto mutation_annotation_test_proto3_optional_new.proto)
	
	# protobuf_generate_cpp(
	# 	PROTO_SRCS PROTO_HDRS
	# 	mutation_annotation_test_proto3_optional_new.proto)

	add_executable(panmat-utils
		${PANMAT_SRCS}
		${PROTO_SRCS}
		${PROTO_HDRS}
		)

endif()

add_dependencies(panmat-utils libminimap2.a)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/clean-minimap)

TARGET_COMPILE_OPTIONS(panmat-utils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
TARGET_LINK_LIBRARIES(panmat-utils PRIVATE stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z)   #${Boost_LIBRARIES}  ZLIB::ZLIB) # OpenMP::OpenMP_CXX)
target_include_directories(panmat-utils PUBLIC "${PROJECT_BINARY_DIR}")