file(GLOB SOURCES "*.cpp")
add_library(compile_tests ${SOURCES})
target_compile_features(compile_tests PRIVATE cxx_std_17)
target_include_directories(compile_tests PUBLIC ${CMAKE_SOURCE_DIR}/include)
