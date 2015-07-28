add_custom_target ( ballet_test )

add_executable( auto/corelib/kernel/test_object
auto/corelib/kernel/test_object.cpp )
target_link_libraries ( auto/corelib/kernel/test_object ballet armadillo )
add_dependencies( ballet_test auto/corelib/kernel/test_object )

add_executable( auto/corelib/tools/test_logical
auto/corelib/tools/test_logical.cpp )
target_link_libraries ( auto/corelib/tools/test_logical ballet armadillo )
add_dependencies( ballet_test auto/corelib/tools/test_logical )

add_executable( auto/corelib/tools/test_random
auto/corelib/tools/test_random.cpp )
target_link_libraries ( auto/corelib/tools/test_random ballet armadillo )
add_dependencies( ballet_test auto/corelib/tools/test_random )

add_executable( auto/comm/tools/test_bi2de
auto/comm/tools/test_bi2de.cpp )
target_link_libraries ( auto/comm/tools/test_bi2de ballet armadillo )
add_dependencies( ballet_test auto/comm/tools/test_bi2de )
