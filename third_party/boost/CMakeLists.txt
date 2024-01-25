file(GLOB_RECURSE bundled_boost_src_files ${PROJECT_SOURCE_DIR}/third_party/boost/libs *.cpp)
add_compile_options(-Wno-effc++ -Wno-shadow)
add_library(bundled_boost ${bundled_boost_src_files})

