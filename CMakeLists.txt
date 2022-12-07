cmake_minimum_required(VERSION 3.15)
project(csr Fortran)

set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(csrmodule STATIC )
set(srcdir src)

option(MKLUSE "Use MKL or not" OFF)
if(MKLUSE)
set(MKL_INTERFACE lp64)
set(MKL_INTERFACE_LAYER "_lp64")
find_package(MKL REQUIRED)
set(lapacklink MKL::MKL)
target_compile_definitions(csrmodule PUBLIC MKL )
else()
endif()

set(linklibraries ${lapacklink})


target_sources(csrmodule 
PRIVATE
${srcdir}/csrmodule.F90
)

target_link_libraries(csrmodule PRIVATE ${linklibraries} )


add_executable(test)

target_sources(test 
PRIVATE
${srcdir}/main.f90
)

target_link_libraries(test csrmodule )




