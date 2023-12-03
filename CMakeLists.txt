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

# armadillo
include_directories(D:/lib/armadillo-12.6.4/include)

# LAPACK
include_directories(D:/lib/lapack/include)
link_directories(D:/lib/lapack/lib)
file(GLOB LAPACK "D:/lib/lapack/lib/*.dll")

link_directories(D:/lib/OpenBLAS-0.3.24-x64/lib)
include_directories(D:/lib/OpenBLAS-0.3.24-x64/include)

add_executable(source
	${SOURCES}
	${HEADS}
	${UIS}
)

target_link_libraries(source D:/lib/armadillo-12.6.4/build/libarmadillo.dll)
target_link_libraries(source D:/lib/OpenBLAS-0.3.24-x64/bin/libopenblas.dll)
target_link_libraries(source ${LAPACK})
target_link_libraries(source Qt5::Widgets)
