cmake_minimum_required(VERSION 3.0.0)
project(PVPSVMC VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(PVPSVMC main.cpp DJ.cpp helper.cpp PVHSS.cpp SVM.cpp)
target_link_libraries(PVPSVMC ntl m gmp pthread relic)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)