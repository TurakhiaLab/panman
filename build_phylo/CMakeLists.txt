cmake_minimum_required (VERSION 3.8) 

project(core_likelihood)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -O3 -pipe -g -ggdb3 -Wno-unused-function  -Wno-deprecated-declarations -lboost_filesystem -lboost_system -lboost_program_options") 
set(CMAKE_INCLUDE_CURRENT_DIR ON) 


# FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


include(${TBB_DIR}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)



add_executable(core_likelihood_fp64
			  ../likelihood_src/llk_fp64.cpp
			  )


add_executable(core_likelihood_fp32
			  ../likelihood_src/llk_fp32.cpp
			  )


TARGET_LINK_LIBRARIES(core_likelihood_fp64 PRIVATE ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS}) #${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
TARGET_LINK_LIBRARIES(core_likelihood_fp32 PRIVATE ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS}) #${Boost_LIBRARIES}  ) # OpenMP::OpenMP_CXX)
