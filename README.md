# CS 112 Final Project

Before compiling, you must change the bunny_path in main.cpp to the absolute path of your own directory where the bunny10k.ply file exists. A example is there.

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

Or, for window users make does not exist so visual studio is needed. To compile for windows, go into the build directory, find a file name example.sln and open it with visual studio. Change the build to release then build solution.

## Run

From within the `build` directory just issue:

    ./example

For window users, the result should be an executable "example.exe" in the release folder in the build directory.

A glfw app should launch displaying a 10k triangles bunny.
