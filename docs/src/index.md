# DCAUtils.jl documentation

```@meta
CurrentModule = DCAUtils
```

This package contains some utilities for [Direct Coupling Analysis](https://en.wikipedia.org/wiki/Direct_coupling_analysis),
an unsupervised technique to analyse Multiple Sequence Alignments of
protein families.

The code is written in [Julia](http://julialang.org). It requires Julia `1.5` or later.

## Installation and Usage

To install the package, use Julia's package manager: from the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
(v1.5) pkg> add DCAUtils
```

Then load it with:

```
julia> using DCAUtils
```

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

