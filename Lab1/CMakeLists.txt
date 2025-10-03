cmake_minimum_required(VERSION 3.10)
project(SplineLab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Исходные файлы
set(SOURCES
    src/main.cpp
    src/Point.cpp
    src/Grid_Generator.cpp
    src/Cubic_Spline.cpp
    src/Numerical_Differentiator.cpp 
)

# Исполняемый файл
add_executable(SplineLab ${SOURCES})

# Настройки компиляции
target_include_directories(SplineLab PRIVATE src)