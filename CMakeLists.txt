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

function(custom_executable directory target)
  add_executable(${target} ${directory}/${target}.cpp)
  target_link_libraries(${target} ${PROJECT_NAME})
  install(
    TARGETS ${target}
    RUNTIME DESTINATION bin/${PROJECT_NAME}
  )
endfunction()

custom_executable(examples access_element)

find_package(GTest QUIET)
if (GTest_FOUND)
    custom_executable(test test_matrix)
    target_link_libraries(test_matrix ${GTEST_LIBRARIES} pthread)
    target_include_directories(test_matrix PUBLIC ${GTEST_INCLUDE_DIRS})
    add_test(NAME test_matrix COMMAND test_matrix)
else (GTest_FOUND)
    MESSAGE("GTest NOT FOUND, tests were not built.")
endif (GTest_FOUND)
