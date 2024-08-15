# Exploiting New Properties of String Net Frequency for Efficient Computation

[![Static Badge](https://img.shields.io/badge/Conference-CPM_2024-blue
)](https://cpm2024.github.io)
[![Static Badge](https://img.shields.io/badge/DOI-10.4230/LIPICS.CPM.2024.16-yellow
)](https://doi.org/10.4230/LIPIcs.CPM.2024.16)

## Introduction

This repository contains the code for the experiments in the following paper:

> Peaker Guo, Patrick Eades, Anthony Wirth, and Justin Zobel: *Exploiting new
properties of string net frequency for efficient computation*. In: 35th Annual Symposium on Combinatorial Pattern Matching (CPM 2024).


## Requirements
* `CMake 3.24+`
* `GCC` or `Clang`


## Building and running

Create the build directory:
```sh
$ mkdir build 
$ cd build
```

Detect the C++ compiler and create build files:
```sh
$ cmake ..
```

Run the build:
```sh
$ cmake --build .
```

Run the program:
```sh
$ ./main
```
