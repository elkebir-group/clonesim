cmake_minimum_required( VERSION 2.8.12 )
project(simcancerphylo)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set(visualizephylo_src
  src/arg_parser.cpp
  src/visualizephylomain.cpp
  src/basetree.cpp
  src/cnatree.cpp
  src/genotypetree.cpp
  src/utils.cpp
  src/phylogeny.cpp
)

set(visualizephylo_hdr
    src/basetree.h
    src/cnatree.h
    src/genotypetree.h
    src/utils.h
    src/phylogeny.h
)

set( mergegenotypetrees_src
    src/mergegenotypetreesmain.cpp
  src/basetree.cpp
  src/genotypetree.cpp
  src/genotypegraph.cpp
  src/cnatree.cpp
  src/cnagraph.cpp
  src/utils.cpp
  src/basematrix.cpp
)

set( mergegenotypetrees_hdr
  src/genotypetree.h
  src/genotypegraph.h
  src/basetree.h
  src/cnatree.h
  src/cnagraph.h
  src/utils.h
  src/basematrix.h
)

set( generategenotypetrees_src
  src/arg_parser.cpp
  src/generategenotypetreesmain.cpp
  src/cnatree.cpp
  src/cnagraph.cpp
  src/basetree.cpp
  src/genotypetree.cpp
  src/genotypegraph.cpp
  src/utils.cpp
  src/basematrix.cpp
)

set( generategenotypetrees_hdr
  src/arg_parser.h
  src/basetree.h
  src/cnatree.h
  src/cnagraph.h
  src/genotypetree.h
  src/genotypegraph.h
  src/utils.h
  src/basematrix.h
)

set( generatecnatrees_src
    src/arg_parser.cpp
    src/generatecnatreesmain.cpp
    src/cnatree.cpp
    src/cnagraph.cpp
    src/genotypetree.cpp
    src/utils.cpp
    src/basematrix.cpp
    src/basetree.cpp
    )

set( generatecnatrees_hdr
    src/arg_parser.h
    src/cnatree.h
    src/cnagraph.h
    src/genotypetree.h
    src/utils.h
    src/basematrix.h
    src/basetree.h
    )

set( simulate_src
  src/arg_parser.cpp
  src/simulatemain.cpp
  src/basetree.cpp
  src/cnatree.cpp
  src/genotypetree.cpp
  src/genotypegraph.cpp
  src/cnagraph.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/phylogeny.cpp
  src/singlecellgeneration.cpp
)

set( simulate_hdr
  src/arg_parser.h
  src/cnatree.h
  src/basetree.h
  src/genotypetree.h
  src/genotypegraph.h
  src/cnagraph.h
  src/utils.h
  src/basematrix.h
  src/phylogeny.h
  src/singlecellgeneration.h
        )

set ( generatesinglecells_src
        src/arg_parser.cpp
        src/singlecellgeneration.cpp
        src/utils.cpp
        src/generatesinglecellsmain.cpp
        src/generatesinglecellsmain.cpp
)

set ( generatesinglecells_hdr
        src/arg_parser.h
        src/utils.h
        src/singlecellgeneration.h
        src/generatesinglecellsmain.cpp
)

set ( pharming_src
        src/arg_parser.cpp
        src/basetree.cpp
        src/utils.cpp
        src/clonaltree.cpp
        src/genotypetree.cpp
        src/cnatree.cpp)

set ( pharming_hdr
        src/arg_parser.h
        src/utils.h
        src/clonaltree.h
        src/genotypetree.h
        src/cnatree.h)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set( Boost_USE_STATIC_LIBS ON )

find_package( Boost COMPONENTS random )

MESSAGE( STATUS "found boost:"
  "\nINCLUDE: ${Boost_INCLUDE_DIRS}"
  "\nLIB: ${Boost_LIBRARIES}"
)

include_directories( "${LIBLEMON_ROOT}/include" ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} "${LIBLEMON_ROOT}/lib" ${GUROBI_LIBRARIES} ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} )

set( CommonLibs
  m
  pthread
  ${Boost_LIBRARIES}
)

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( mergegenotypetrees ${mergegenotypetrees_src} ${mergegenotypetrees_hdr} )
target_link_libraries( mergegenotypetrees ${CommonLibs} )

add_executable( generategenotypetrees ${generategenotypetrees_src} ${generategenotypetrees_hdr} )
target_link_libraries( generategenotypetrees ${CommonLibs} )

add_executable( generatecnatrees ${generatecnatrees_src} ${generatecnatrees_hdr} )
target_link_libraries( generatecnatrees ${CommonLibs} )

add_executable( simulate ${simulate_src} ${simulate_hdr} )
target_link_libraries( simulate ${CommonLibs} )

add_executable( visualizephylo ${visualizephylo_src} ${visualizephylo_hdr} )
target_link_libraries( visualizephylo ${CommonLibs} )

add_executable( generatesinglecells ${generatesinglecells_src} ${generatesinglecells_hdr}
        src/generatesinglecellsmain.cpp)
target_link_libraries( generatesinglecells ${CommonLibs})

add_executable( pharming ${pharming_src} ${pharming_hdr}
        src/clonaltree.cpp)
target_link_libraries( pharming ${CommonLibs})