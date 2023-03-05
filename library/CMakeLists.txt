set(TARGET_NAME celib)

# FIND OPENMP
find_package(OpenMP)

find_package(Ceres)

find_package(yaml-cpp)

find_package(OpenCV REQUIRED)

find_package(Boost COMPONENTS iostreams REQUIRED)

file(GLOB source src/*/*.cpp)

add_library( ${TARGET_NAME} STATIC ${source} )
set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_compile_definitions( ${TARGET_NAME} PUBLIC ${PCL_DEFINITIONS})
target_include_directories( ${TARGET_NAME} PUBLIC
    include
    ${CERES_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/external_library)

target_link_libraries( ${TARGET_NAME}
    ${CERES_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
    ${catkin_LIBRARIES}
    ${OpenMP_CXX_LIBRARIES}
    ${OpenCV_LIBS}
    libgnuplot-iostream
    cilantro
    )

target_compile_options( ${TARGET_NAME} PRIVATE
    ${OpenMP_CXX_FLAGS})
