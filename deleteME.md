## Installing / Getting started
```shell
git clone https://github.com/vzendejasl/CS179.git
```

## Installation Requirements
The following support modules (or newer) are required to run the code,

1. gcc version 5.0.0 
2. cmake version 3.9.1 
3. nvcc 9.1 (Nvidia toolkit)

The Nvidia toolkit can be downloaded from https://developer.nvidia.com/cuda-toolkit .
You can check your cmake version and gcc version with

```shell
gcc --version
vcc --version
```
## Building
To compile, in the `232aCFDCode_CUDA` directory, simply run
```shell
mkdir build
cd build
```

If you are degguing/adding new features to the code then,
```shell
cmake -DCMAKE_BUILD_TYPE=Debug ../
```
otherwise,
```shell
cmake -DCMAKE_BUILD_TYPE=Release ../
```
and finally
```shell
cmake --build .
```

## Running
Once the makefile is created `build`,

```shell
make 
```
and you will see the executable `2dSolverNS`. 

To run the code run,
```shell
./2dSolverNS
```
You should see the following output in your terminal,

```shell
Running CPU Code...
Iteration # 0
Iteration # 1
Iteration # 2
Iteration # 3
Iteration # 4
Time taken by CPU: 20680 milliseconds
---------------------------------
Running GPU Code...
Iteration # 0
Iteration # 1
Iteration # 2
Iteration # 3
Iteration # 4
Time taken by GPU: 788 milliseconds
```

You can control the number of time steps taken by the code by modifying `int num_iter = 5;` with a larger value in `main.cpp`.

## Visualization
You should see two text files in the same directory as `2dSolverNS`, `CPU_output.txt` and `GPU_output.txt`. These text files are the CPU and GPU solutions at the final time step. In the main directoy you will find a MATLAB script `PlotOutput.m`. You can run the scipt as any other MATLAB program with the exeption that the number of grid points need to be specified and need to match the number of grid points used in the calculations. The example shown below, `NX = NY = 2048`. 



Usage Instructions - Victor 
- What installation steps are necessary (if any)?
- How do we run the program and see output? (please include a few demo scripts!)

Project Description - Alex
- What does the program do?

Results - Victor
- What should we expect to see?
Please include some screenshots!

Performance Analysis - Alex
- How much better is the GPU version? If it's not better, why not? (This sometimes happens, oops, when there is a big unexpected conflict on the GPU. Many times this can be resolved, even with a small change, but perhaps an issue like this might not be findable in a four week project!)
- Are there things that could be improved?
