cmake_minimum_required (VERSION 3.26) 

project(panmat-utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17 -Wall -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 

find_package(Protobuf REQUIRED CONFIG) 
include_directories(${Protobuf_INCLUDE_DIRS})
if(Protobuf_FOUND)
	message(STATUS "Found Protocol Buffers: ${Protobuf_VERSION}")
endif()

FIND_PACKAGE(Boost COMPONENTS system filesystem program_options unit_test_framework iostreams date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS arch=arm64 tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

file(GLOB PANMAT_SRCS "src/main.cpp" "src/PangenomeMAT.cpp" "src/PangenomeMAT.hpp")
file(GLOB PANMAT_LIB_SRCS "src/PangenomeMAT.cpp" "src/PangenomeMAT.hpp")

# targets
add_executable(panmat-utils ${PANMAT_SRCS})

protobuf_generate(LANGUAGE cpp 
	TARGET panmat-utils 
	PROTOS mutation_annotation_test_proto3_optional_new.proto)

target_compile_options(panmat-utils PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_libraries(panmat-utils PRIVATE ${JSONCPP_LINK_LIBRARIES} ${spoa_DIR}/build/lib/libspoa.a stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} protobuf::libprotobuf z)
target_include_directories(panmat-utils PUBLIC "${PROJECT_BINARY_DIR}" ${spoa_DIR}/include)



add_library(libpanmat SHARED ${PANMAT_LIB_SRCS})
protobuf_generate(LANGUAGE cpp 
	TARGET libpanmat
	PROTOS mutation_annotation_test_proto3_optional_new.proto)
target_link_libraries(libpanmat PRIVATE ${JSONCPP_LINK_LIBRARIES} ${spoa_DIR}/build/lib/libspoa.a stdc++ ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} protobuf::libprotobuf z)
target_include_directories(libpanmat PUBLIC "${PROJECT_BINARY_DIR}" ${spoa_DIR}/include)
set_target_properties(libpanmat PROPERTIES PUBLIC_HEADER src/PangenomeMAT.hpp)
export(TARGETS libpanmat FILE libpanmatConfig.cmake)

configure_file(src/PangenomeMAT.hpp include/PangenomeMAT.hpp COPYONLY)
install(TARGETS libpanmat
	ARCHIVE DESTINATION "lib"
	LIBRARY DESTINATION "lib"
	PUBLIC_HEADER DESTINATION "include"
	)
