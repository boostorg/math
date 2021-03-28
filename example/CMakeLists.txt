file(GLOB SOURCES "*.cpp")
add_library(examples ${SOURCES})
target_compile_features(examples PRIVATE cxx_std_17)
target_include_directories(compile_tests PUBLIC ${CMAKE_SOURCE_DIR}/include)
