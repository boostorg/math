file(GLOB SOURCES "*.cpp")
link_directories(${CMAKE_SOURCE_DIR}/include)
add_library(examples ${SOURCES})
target_compile_features(examples PRIVATE cxx_std_17)
