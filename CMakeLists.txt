# cmake版本要求
cmake_minimum_required(VERSION 3.0)

# 项目名称
project(generatorEvent)

# 设置编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# 找root
find_package(ROOT REQUIRED)

# 设置头文件路径
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

# 设置库文件路径
link_directories(${ROOT_LIBRARY_DIR})
link_directories(${PROJECT_SOURCE_DIR}/lib)

# 设置源文件
set(SOURCE_FILES src/*.cpp)

# 生成可执行文件
add_executable(generatorEvent generatorEvent.cpp ${SOURCE_FILES})

# 链接库文件
target_link_libraries(generatorEvent ${ROOT_LIBRARIES}  )