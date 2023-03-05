# FIND BOOST
find_package(Boost COMPONENTS program_options REQUIRED)

find_package(yaml-cpp)


# IN2LAAMA demonstration
add_executable(gp_comet src/egp_demo.cpp)
set_target_properties(gp_comet PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories( gp_comet PUBLIC
    ${PROJECT_SOURCE_DIR}/external_library/
    ${PROJECT_SOURCE_DIR}/library/
    ${Boost_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    )

target_link_libraries( gp_comet PUBLIC
    celib
    ${Boost_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
    )
