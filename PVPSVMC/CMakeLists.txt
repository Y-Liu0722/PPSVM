cmake_minimum_required(VERSION 3.0.0)
project(PVPSVMC VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(PVPSVMC
        main.cpp
        helper.cpp
        helper.h
        SVM.cpp
        SVM.h
        PKE.cpp
        PKE.h
        BKS.cpp
        BKS.h
        bench.cpp
        bench.h)

target_link_libraries(PVPSVMC ntl m gmp pthread relic helib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)