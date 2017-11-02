# OpenLSTO
**Open Source Level Set Topology Optimization Software**

OpenLSTO is an open-source software for level set based structural topology optimization, written in C++. Developed and Maintained by researchers at UC San Diego and Cardiff University.

The topology optimization method deals with the problem of determining the optimal layout of material and connectivity inside a design domain in order to obtain the best structural performance. Level Set Method was originally developed as a mathematical tool for tracking the motion of interfaces. Its natural handling of topological changes coupled with a clear and smooth interface representation led to its use for structural topology optimization.

This first light version implements the M2DO lab’s level set topology optimization method to solve the problem of minimizing compliance (maximizing stiffness) under a volume constraint in a 2D design domain. It has been designed with ease of installation and use in mind. This means that, wherever possible, a conscious effort was made to develop in-house code components rather than relying on third-party packages or libraries. In simple cases (serial version with no external libraries), the code can be compiled and executed with just a C++ compiler.

Extended capabilities using externally provided software are to be included in the next, more advanced version.

## Installation
-  Download

To clone the repository to the local machine,(git client) 
```git clone https://github.com/M2DOLab/OpenLSTO.git```

The Zip file is also avaiable with the following static link: 
[Download link](http://m2do.ucsd.edu/static/zip/OpenLSTO-v0.1.zip)

- CLI environment

In general, all OpenLSTO execution occurs via command line arguments within a terminal. For Unix/Linux or Mac OS X users, the native terminal applications are needed. 

- Execution

Compile & Build: 
> make 

Execution
> ./a.out

Clean existing compilation files:
> make clean


- External dependency

Users of OpenLSTO need a data visualization tool to post-process solution files. The software currently supports .vtk output format natively read by [*ParaView*](https://www.paraview.org/).

## Licensing
OpenLSTO is available for download under the [Apache V. 2.0](https://www.apache.org/licenses/LICENSE-2.0) license. Please refer to the License page for terms and conditions.

## Contributors

Prof. Hyunsun Alicia Kim (UCSD, Cardiff U.)

Dr. Sandilya Kambampati (UCSD)

Dr. Lester Hedges (Lawrence Berkeley National Lab)

Dr. Zongliang Du (UCSD)

Dr. Renato Picelli (Cardiff U.)

Dr. Scott Townsend (Cardiff U.)

Dr. Xiao-Yi Zhou (Cardiff U.)

Dr. Hayoung Chung (UCSD)

Ms. Carolina Jauregui (UCSD)
