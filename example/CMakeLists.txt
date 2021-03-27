file(GLOB SOURCES "*.cpp")
add_library(examples ${SOURCES})
target_compile_features(examples PRIVATE cxx_std_17)
