cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

if(NOT (UNIX))
    message(FATAL_ERROR "This project is intended for UNIX systems only.")
endif()

project(uwb)

include(ExternalProject)

execute_process(COMMAND rm -rf ../wavelib)
execute_process(COMMAND rm -rf ../iir1)

execute_process(COMMAND git clone "https://github.com/rafat/wavelib" ${CMAKE_CURRENT_SOURCE_DIR}/wavelib)
execute_process(COMMAND git clone "https://github.com/berndporr/iir1" ${CMAKE_CURRENT_SOURCE_DIR}/iir1)

add_subdirectory(wavelib)
add_subdirectory(iir1)

add_executable(test_write test_write.cpp uwb.h uwb.cpp)
add_executable(test_read test_read.cpp uwb.cpp uwb.h)
add_executable(test_demodulate test_demodulate.cpp uwb.cpp uwb.h)
add_executable(test_background_subtraction test_background_subtraction.cpp uwb.cpp uwb.h)
add_executable(test_distance test_distance.cpp uwb.h uwb.cpp distance.h distance.cpp)
add_executable(test_presence test_presence.cpp uwb.h uwb.cpp)
add_executable(test_decimation test_decimation.cpp uwb.h uwb.cpp)

target_link_libraries(test_distance wavelib iir::iir)
target_link_libraries(test_presence wavelib iir::iir)
target_link_libraries(test_write wavelib iir::iir)
target_link_libraries(test_read wavelib iir::iir)
target_link_libraries(test_demodulate wavelib iir::iir)
target_link_libraries(test_background_subtraction wavelib iir::iir)
target_link_libraries(test_decimation wavelib iir::iir)

set_target_properties(test_write test_read test_demodulate test_background_subtraction test_distance test_presence
   PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
)


