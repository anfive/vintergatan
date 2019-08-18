# Vintergatan CMake file.
# Author: Andrea Ferrario

cmake_minimum_required(VERSION 3.10)

project(vintergatan VERSION 1.0.0)

add_executable(vintergatan
	settings.h
	graphics.cpp graphics.h
	main.cpp
	octree.cpp octree.h
	octree_tests.cpp
	simulation.cpp simulation.h
	testcases.cpp testcases.h
)
set_property(TARGET vintergatan PROPERTY CXX_STANDARD 17)

# Manage dependencies.

# OpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
	message(FATAL_ERROR "OpenMP not found")    
endif()
target_link_libraries(vintergatan PUBLIC OpenMP::OpenMP_CXX)


# OpenGL
find_package(OpenGL REQUIRED)
if(NOT OpenGL_FOUND)
	message(FATAL_ERROR "OpenGL not found")
endif()
target_include_directories(vintergatan PUBLIC ${OPENGL_INCLUDE_DIR})
target_link_libraries(vintergatan PUBLIC ${OPENGL_gl_LIBRARY})

# GLEW
# Link to the static library so we don't have to manage the dll on windows.
set(glew-cmake_BUILD_SHARED OFF)
set(glew-cmake_BUILD_STATIC ON)
set(glew-cmake_BUILD_SINGLE_CONTEXT ON)
set(glew-cmake_BUILD_MULTI_CONTEXT OFF)
add_subdirectory(glew-cmake)
target_compile_definitions(libglew_static PUBLIC GLEW_STATIC) # A bug in glew-cmake? It does not define GLEW_STATIC for the static library target.
target_include_directories(vintergatan PRIVATE libglew_static)
target_link_libraries(vintergatan PUBLIC libglew_static)

# GLFW
set(GLFW_BUILD_DOCS False)
set(GLFW_INSTALL False)
add_subdirectory(glfw)
target_link_libraries(vintergatan PUBLIC glfw)