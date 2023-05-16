cmake_minimum_required(VERSION 3.0)
project(Lift-Generation-Model)

# Set the C standard to use (optional)
set(CMAKE_C_STANDARD 11)

# Add your source files
set(SOURCES src/lift_gen.c)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

