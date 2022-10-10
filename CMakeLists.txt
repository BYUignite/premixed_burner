cmake_minimum_required(VERSION 3.15)

project(cantera_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(burner.x)
target_sources(burner.x PRIVATE src/burner.cc)

find_library(SOOT sootModel HINTS "/Users/dol4/codes/sootlib/lib")
find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)
if(APPLE)
    find_library(ACCELERATE Accelerate)
endif()

target_include_directories(burner.x PRIVATE "/usr/local/include" 
                                            "/opt/homebrew/include"
                                            "/Users/dol4/codes/sootlib/include/sootlib")
target_link_libraries(burner.x ${CANTERA} ${CVODES} ${SOOT} ${ACCELERATE})

install(TARGETS burner.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
