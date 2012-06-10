
This repository hosts the code of **Symmetrix**.

# What is Symmetrix?

**Symmetrix** -abreviated **SYMTRX**- is an open-source library that provides specific algorithms and benchmark for exotic matrices.

**Symmetrix** do not aim at challenging large projects like [BLAS](http://www.netlib.org/blas/)/[LAPACK](http://www.netlib.org/lapack/), [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) or even [MKL](http://software.intel.com/en-us/articles/intel-mkl/). Its purpose is simply to provide a collection of efficient algorithms for improving the speed and saving resources when dealing with exotic matrices. For this reason all vectors and matrices are assumed to be dense, and the related algorithms are rather naive. You should keep using the previously-mentioned libraries for the usual vector-matrix operations. For exotic matrices though, all algorithms tend to be optimised and benchmarked so that they can be used in large problems.

# Current features

* C library
	* Vectors -
	* Square matrices - product, trace, traceproduct
	* Centrosymmetric matrices - product, trace, traceproduct
	* Bisymmetric matrices -

# Remarks

The documentation can be generated with doxygen.

We are (I am) very open to remarks, contributions and feedback!

Visit [Boris Leistedt's blog](http://ixkael.com/blog) for maths-related documentation, announcements and applications of the library in Astrophysics.


