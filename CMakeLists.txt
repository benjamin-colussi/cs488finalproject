cmake_minimum_required(VERSION 3.31)
project(finalproject)

set(CMAKE_CXX_STANDARD 20)

add_executable(finalproject main.cpp)



#find_package(OpenMP REQUIRED)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(finalproject PRIVATE OpenMP::OpenMP_CXX)
#endif()



#==> Caveats
#libomp is keg-only, which means it was not symlinked into /usr/local,
#because it can override GCC headers and result in broken builds.
#
#For compilers to find libomp you may need to set:
#export LDFLAGS="-L/usr/local/opt/libomp/lib"
#export CPPFLAGS="-I/usr/local/opt/libomp/include"
#==> Summary
#🍺  /usr/local/Cellar/libomp/20.1.7: 9 files, 1.7MB