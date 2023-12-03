# windeployqt source.exe
# windeployqt math software.exe
cmake_minimum_required(VERSION 3.5)

#add_definitions(-DDEBUG)
#add_definitions(-DPRINT_LATEX)
#add_definitions(-DZSC_TIME)
#add_definitions(-DZSC_ERROR)
#add_definitions(-DDEBUG41)
#add_definitions(-DTIME41)
#add_definitions(-DDEBUG42)
#add_definitions(-DDEBUG43)
#add_definitions(-DDEBUG43E1)
#add_definitions(-DDEBUG44)
#add_definitions(-DDEBUG51)
#add_definitions(-DDEBUG52)
#add_definitions(-DTIME52)
#add_definitions(-DDEBUG53)
#add_definitions(-DDEBUG61)
#add_definitions(-DDEBUG611)
#add_definitions(-DDEBUG28)
#add_definitions(-DDEBUG281)
#add_definitions(-DDEBUG282)
#add_definitions(-DDEBUG62)
#add_definitions(-DDEBUG331)
#add_definitions(-DDEBUG31)
#add_definitions(-DDEBUG34)
add_definitions(-DDEBUG32)

project(source LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLAPACK_DLL")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLAPACK_DLL")

find_package(Qt5 COMPONENTS Widgets REQUIRED)

file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE UIS "src/*.ui")
file(GLOB_RECURSE HEADS "include/*.h")
include_directories(${PROJECT_SOURCE_DIR}/include/Qt)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex1)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex2)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex3)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex4)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex5)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex6)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex7)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex8)
include_directories(${PROJECT_SOURCE_DIR}/include/Ex9)
include_directories(${PROJECT_SOURCE_DIR}/include)

# 以下请根据自己的安装路径进行修改
#set(LIBS
#	# armadillo库所需的静态库
#	liblapack.lib
#	liblapacke.lib
#	libblas.lib
#)
#set(DLLS
#	D:/lib/lapack/liblapack.dll
#	D:/lib/lapack/liblapacke.dll
#	D:/lib/lapack/libblas.dll
#)

# armadillo
#add_subdirectory(D:/lib/armadillo-12.6.4 ./bin)
include_directories(D:/lib/armadillo-12.6.4/include)
#file(GLOB_RECURSE ARMADILLO_SOURCES "D:/lib/armadillo-12.6.4/*.cpp")
#find_package(ARMADLIIO_LIBRARY CONFIG REQUIRED NAMES Armadillo HINTS D:/lib/armadillo-12.6.4/build)
#find_package(armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})
#link_directories(D:/lib/armadillo-12.6.4/build)
#link_libraries(Armadillo)
#include_directories(${ARMADLIIO_LIBRARY_INCLUDE_DIR})
#link_directories(D:/lib/armadillo-12.6.4/examples/lib_win64)
#file(GLOB LIBS D:/lib/armadillo-12.6.4/examples/lib_win64)
#link_libraries(${LIBS})

# LAPACK
#FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack HINTS D:/lib/armadillo-12.6.4/lib)
include_directories(D:/lib/lapack/include)
link_directories(D:/lib/lapack/lib)
file(GLOB LAPACK "D:/lib/lapack/lib/*.dll")
#link_libraries(${LIBS})
#find_package(LAPACK_LIBRARY REQUIRED NAMES LAPACK HINTS D:/lib/armadillo-12.6.4/cmake_aux/Modules)
#include_directories(${LAPACK_LIBRARY_INCLUDE_DIR})
#target_link_libraries(source PUBILC ${LAPACK_LIBRARY_LIBRARY})

# blas
#find_package(BLAS_LIBRARY REQUIRED NAMES BLAS HINTS D:/lib/armadillo-12.6.4/cmake_aux/Modules)
#include_directories(${BLAS_LIBRARY_INCLUDE_DIR})
#target_link_libraries(source PUBILC ${BLAS_LIBRARY_LIBRARY})
link_directories(D:/lib/OpenBLAS-0.3.24-x64/lib)
include_directories(D:/lib/OpenBLAS-0.3.24-x64/include)

# MathJax
#set(MATHJAX_DIR "D:/lib/MathJax")
#find_path(MathJax_INCLUDE_DIRS
#  NAMES MathJax.js
#  HINTS ${MATHJAX_DIR}
#)
#add_library(MathJax INTERFACE)
#target_include_directories(MathJax INTERFACE ${MathJax_INCLUDE_DIRS})
#target_link_libraries(source PRIVATE MathJax)

# OpenCV
#include_directories("D:/lib/opencv/build/include")
#file(GLOB OPENCV_LIB "D:/lib/opencv/build/x64/vc16/lib/*.lib")
#link_libraries(${OPENCV_LIB})

#target_link_libraries(source PUBLIC ${DLLS})
#link_libraries(${LIBS})
#target_link_libraries(source PUBLIC ${LIBS})

add_executable(source
	${SOURCES}
	${HEADS}
	${UIS}
)

#target_link_libraries(source ${ARMADLIIO_LIBRARY_LIBRARY})
#target_link_libraries(source PUBLIC ${LIBS})
#link_directories(${DLLS})

#add_definitions(-DLAPACK_DLL)  # 添加编译选项

#set(LAPACK_LIB_PATH D:/lib/lapack/lib)  # 设置链接库的路径

# 手动添加链接库
#if(WIN32)
#	# Windows平台
#	foreach(LAPACK_LIB "lapack" "blas")
#		if(MSVC)
#			# 使用Visual Studio编译器
#			target_link_libraries(${PROJECT_NAME} ${LAPACK_LIB_PATH}/${LAPACK_LIB}.lib)
#		else()
#			# 使用其他编译器
#			target_link_libraries(${PROJECT_NAME} ${LAPACK_LIB_PATH}/${LAPACK_LIB}.dll)
#		endif()
#	endforeach()
#elseif(APPLE)
#	# macOS平台
#	foreach(LAPACK_LIB "lapack" "blas")
#		target_link_libraries(${PROJECT_NAME} ${LAPACK_LIB_PATH}/${LAPACK_LIB}.a)
#	endforeach()
#else()
#	# 其他Linux平台
#	foreach(LAPACK_LIB "lapack" "blas")
#		target_link_libraries(${PROJECT_NAME} ${LAPACK_LIB_PATH}/${LAPACK_LIB}.a)
#	endforeach()
#endif()

#target_link_libraries(source PUBLIC ${ARMADILLO_LIB_DIRS} armadillo ${LAPACK_LIBRARY})
#target_link_libraries(source PUBLIC ${ARMADILLO_LIB_DIRS} armadillo -O1 -llapack -lblas)
#target_link_libraries(source PUBLIC ${ARMADILLO_LIB_DIRS} armadillo)
#target_link_libraries(source PUBLIC armadillo)
target_link_libraries(source D:/lib/armadillo-12.6.4/build/libarmadillo.dll)
target_link_libraries(source D:/lib/OpenBLAS-0.3.24-x64/bin/libopenblas.dll)
#target_link_libraries(source ${LIBS})
target_link_libraries(source ${LAPACK})
target_link_libraries(source Qt5::Widgets)
