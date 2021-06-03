# Performance Analysis in HPC

No matter how complex your scientific application is, with tens of millions of code lines and  humonguous data structures, at the core of it there are always some quite basic operations that are performed in a CPU/GPU. These are, most often, arithmetic operations and array operations. The initial step to write simple and efficient programs in languages that are somewhat close to the underlying machine instructions, such as Fortran and C++, is therefore understanding how the hardware you are using performs these operations and, also, how long it takes to perform them. The idea of this project is presenting some very simple code that will give you a taste of what your processing unit is capable of doing, so you can find your way to write an efficient code from the get-go. 


## basics
In this folder, there are straightforward codes that time arithmetic and array operations. For the former, the focus is on comparing single precision (sp) vs. double precision (dp) calculations. You should compile and run these without having any troubles. Note that different compilers and compilation flags may give you different results. I encourage you to explore these. In general, sp and dp are pretty close together except for when divisions are required. For the latter, the idea is noticing what is the correct order to perform nested loops that enclose tensor operations. This depends on how your compiler uses contiguous memory allocations, which is fundamentally different for Fortran (column-major) and C++ (row-major). Using the correct loop order saves a lot of time, as it dramatically increases cache efficiency.




## CITE-ME
As all other projects in this GitHub, the programs are protected by a license/copyright that requires proper acknowledgement. Please refer to: https://github.com/babreu-ncsa/cite-me.
