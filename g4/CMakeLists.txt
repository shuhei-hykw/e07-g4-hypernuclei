cmake_minimum_required(VERSION 3.16)
project(emulsim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Physics Tree RIO)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
add_executable(emulsim emulsim.cc ${sources})
target_include_directories(emulsim PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(emulsim ${Geant4_LIBRARIES} ROOT::Physics
                      ROOT::Tree ROOT::RIO)

file(GLOB macros ${PROJECT_SOURCE_DIR}/macros/*.mac)
file(COPY ${macros} DESTINATION ${PROJECT_BINARY_DIR})
