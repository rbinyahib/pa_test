cmake_minimum_required(VERSION 2.8)
project(testing)


SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

#Use this on cori
#SET(CMAKE_C_COMPILER cc)
#SET(CMAKE_CXX_COMPILER CC)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
include_directories(.)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#find_package(VTKm REQUIRED)
#include_directories(${VTKm_INCLUDE_DIRS})


find_package(VTKm REQUIRED
  OPTIONAL_COMPONENTS Rendering
  OPTIONAL_COMPONENTS TBB
  )


# find_package(VTKm REQUIRED QUIET
#   OPTIONAL_COMPONENTS Serial OpenGL Rendering GLUT
#   )

# if(VTKm_OpenGL_FOUND AND VTKm_Rendering_FOUND AND VTKm_GLUT_FOUND)
#   add_executable(Rendering_SERIAL Rendering.cxx)
#   target_include_directories(Rendering_SERIAL PRIVATE ${VTKm_INCLUDE_DIRS})
#   target_link_libraries(Rendering_SERIAL ${VTKm_LIBRARIES})
#   target_compile_options(Rendering_SERIAL PRIVATE ${VTKm_COMPILE_OPTIONS})
# endif()



add_executable(Testing TimingsManager.C testing2.cxx)
#add_executable(Testing testing2.cxx)

link_directories(${VTK_LIBRARY_DIRS})
#find_path(DIY_INCLUDE_DIRS diy/types.h)

#include_directories(${DIY_INCLUDE_DIRS})
#include_directories(${SPDLOG_INCLUDE_DIR})

#target_link_libraries(Testing pcomm)
target_include_directories(Testing PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(Testing ${VTKm_LIBRARIES})
target_compile_options(Testing PRIVATE ${VTKm_COMPILE_OPTIONS})

target_link_libraries(Testing ${libraries} ${VTK_LIBRARY} vtkIOCore vtkIOLegacy vtkImagingCore vtkFiltersFlowPaths)

#if(VTKm_TBB_FOUND)
#  add_executable(testing_TBB testingTBB.cxx)
#  target_include_directories(testing_TBB PRIVATE ${VTKm_INCLUDE_DIRS})
#  target_link_libraries(testing_TBB ${VTKm_LIBRARIES})
#  target_compile_options(testing_TBB PRIVATE ${VTKm_COMPILE_OPTIONS})
#endif()

