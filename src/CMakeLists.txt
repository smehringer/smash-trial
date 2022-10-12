cmake_minimum_required (VERSION 3.8)

# A interface to reuse common properties.
# You can add more external include paths of other projects that are needed for your project.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include)
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" STATIC search.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC "${PROJECT_NAME}_interface")
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC raptor_interface)
target_compile_definitions ("${PROJECT_NAME}_lib" PUBLIC "-DRAPTOR_HIBF_HAS_COUNT=1")

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib_j" STATIC jaqquard_dist.cpp)
target_link_libraries ("${PROJECT_NAME}_lib_j" PUBLIC "${PROJECT_NAME}_interface")
target_link_libraries ("${PROJECT_NAME}_lib_j" PUBLIC raptor_interface)
target_compile_definitions ("${PROJECT_NAME}_lib_j" PUBLIC "-DRAPTOR_HIBF_HAS_COUNT=1")

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib_j")

add_executable (matrixify_mash_output matrixify_mash_output.cpp)
target_link_libraries (matrixify_mash_output PRIVATE "${PROJECT_NAME}_lib")
target_link_libraries (matrixify_mash_output PRIVATE "${PROJECT_NAME}_lib_j")

add_executable (sort_matrix sort_matrix.cpp)
target_link_libraries (sort_matrix PRIVATE "${PROJECT_NAME}_lib")
target_link_libraries (sort_matrix PRIVATE "${PROJECT_NAME}_lib_j")
