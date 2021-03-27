file(GLOB SOURCES "*.cpp")
link_directories(${CMAKE_SOURCE_DIR}/include)
add_library(compile_tests ${SOURCES})
target_compile_features(compile_tests PRIVATE cxx_std_17)
