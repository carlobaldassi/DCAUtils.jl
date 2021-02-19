# DCAUtils.jl documentation

```@meta
CurrentModule = DCAUtils
```

This package contains some utilities for [Direct Coupling Analysis](https://en.wikipedia.org/wiki/Direct_coupling_analysis),
an unsupervised technique to analyse Multiple Sequence Alignments of
protein families.

The code is written in [Julia](http://julialang.org). It requires Julia `1.5` or later.

## Installation and Usage

To install the package, use Julia's package manager: from the Julia REPL, type `]` to enter the Pkg
REPL mode and run:

```
(v1.5) pkg> add DCAUtils
```

Then load it with:

```
julia> using DCAUtils
```

## Overview

The functions in this package are written to maximize performance. Most computationally-heavy
functions can use multiple threads (start julia with the `-t` option or set the `JULIA_NUM_THREADS`
environment variable). In most cases such functions need to perform operations over every pair of
sequences or every pair of resiudes; since such operations are symmetric, only the upper-triangular
part is needed, and a custom parallelization scheme is used for this purpose. During these parallel
parts of the code, BLAS parallelization is disabled, and it is restored at the end.

On top of parallelization, big efficiency gains can be had in certain operations by working with a
compressed representation of the data. The general idea is that the multiple sequence alignment is
first parsed and converted into a `Matrix{Int8}` representation (see
[`read_fasta_alignment`](@ref)). Then, internally, functions such as [`compute_weights`](@ref) and
[`compute_dists`](@ref) further explot the assumption that no value larger than `31` will appear,
thus requiring 5 bits at the most: thus, they pack each sequence of $N$ `Int8` values into $⌈N /
12⌉$ `UInt64` values. Efficient bit-wise operations are then applied to compute the Hamming
distance between sequences, processing 12 entries at a time.


## Reference

```@docs
read_fasta_alignment
```
```@docs
remove_duplicate_sequences
```
```@docs
compute_theta
```
```@docs
compute_weights
```
```@docs
compute_dists
```
```@docs
compute_weighted_frequencies
```
```@docs
add_pseudocount
```
```@docs
compute_DI_gauss
```
```@docs
compute_FN
```
