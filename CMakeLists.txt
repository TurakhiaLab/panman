cmake_minimum_required (VERSION 3.26) 

project(panmat-utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++17") 
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Protobuf
find_package(Protobuf REQUIRED CONFIG)
include_directories(${Protobuf_INCLUDE_DIRS})
if(Protobuf_FOUND)
	message(STATUS "Found Protocol Buffers: ${Protobuf_VERSION}")
endif()

# Boost
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIRS})

# TBB
include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS arch=arm64 tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

# minimap2
add_custom_target(libminimap2.a
   COMMAND make arm_neon=1 aarch64=1 libminimap2.a
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/minimap2_src/
)
add_custom_target(clean-minimap
    COMMAND make clean
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/minimap2_src/
)
file(GLOB PANMAT_SRCS "src/align_from_seeds.c" "src/minimap2_src/libminimap2.a" "src/*.cpp" "src/*.hpp")

# targets
add_executable(panmat-utils ${PANMAT_SRCS})

protobuf_generate(LANGUAGE cpp 
	TARGET panmat-utils 
	PROTOS mutation_annotation_test_proto3_optional.proto mutation_annotation_test_proto3_optional_new.proto)

	
add_dependencies(panmat-utils libminimap2.a)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/clean-minimap)
target_compile_options(panmat-utils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_libraries(panmat-utils PRIVATE ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} protobuf::libprotobuf ${CMAKE_SOURCE_DIR}/src/minimap2_src/libminimap2.a z) 
target_include_directories(panmat-utils PUBLIC "${PROJECT_BINARY_DIR}")