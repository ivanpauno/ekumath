cmake_minimum_required(VERSION 3.5)

project(ekumath)

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_compile_options(-Wall -Wextra -Wpedantic)

include_directories(include)

add_library(${PROJECT_NAME} src/matrix.cpp)

install(
  TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

function(custom_executable target)
  add_executable(${target} examples/${target}.cpp)
  target_link_libraries(${target} ${PROJECT_NAME})
  install(
    TARGETS ${target}
    RUNTIME DESTINATION bin/${PROJECT_NAME}
  )
endfunction()

custom_executable(access_element)

find_package(GTest QUIET)
if (GTest_FOUND)
    MESSAGE("Included GTest")
else (GTest_FOUND)
    MESSAGE("GTest NOT FOUND. Searched it in: ${GTEST_INCLUDE_DIRS}")
endif (GTest_FOUND)

