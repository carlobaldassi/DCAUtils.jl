var documenterSearchIndex = {"docs":
[{"location":"#DCAUtils.jl-documentation","page":"Home","title":"DCAUtils.jl documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = DCAUtils","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package contains some utilities for Direct Coupling Analysis, an unsupervised technique to analyse Multiple Sequence Alignments of protein families.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The code is written in Julia. It requires Julia 1.5 or later.","category":"page"},{"location":"#Installation-and-Usage","page":"Home","title":"Installation and Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package, use Julia's package manager: from the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(v1.5) pkg> add DCAUtils","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then load it with:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using DCAUtils","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The functions in this package are written to maximize performance. Most computationally-heavy functions can use multiple threads (start julia with the -t option or set the JULIA_NUM_THREADS environment variable). In most cases such functions need to perform operations over every pair of sequences or every pair of resiudes; since such operations are symmetric, only the upper-triangular part is needed, and a custom parallelization scheme is used for this purpose. During these parallel parts of the code, BLAS parallelization is disabled, and it is restored at the end.","category":"page"},{"location":"","page":"Home","title":"Home","text":"On top of parallelization, big efficiency gains can be had in certain operations by working with a compressed representation of the data. The general idea is that the multiple sequence alignment is first parsed and converted into a Matrix{Int8} representation (see read_fasta_alignment). Then, internally, functions such as compute_weights and compute_dists further explot the assumption that no value larger than 31 will appear, thus requiring 5 bits at the most: thus, they pack each sequence of N Int8 values into N  12 UInt64 values. Efficient bit-wise operations are then applied to compute the Hamming distance between sequences, processing 12 entries at a time.","category":"page"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"read_fasta_alignment","category":"page"},{"location":"#DCAUtils.ReadFastaAlignment.read_fasta_alignment","page":"Home","title":"DCAUtils.ReadFastaAlignment.read_fasta_alignment","text":"read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real) -> Matrix{Int8}\n\nParses a FASTA file containing a multiple sequence alignment, and returns a matrix of integers that represents one sequence per column.\n\nThe mapping between the aminoacid symbols and the integers uses this table:\n\n  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y\n  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20\n\nAny unrecognized capital letter and the gap symbol - are mapped to the value 21; any other symbol or lowercase letter is ignored.\n\nIf a sequence contains a fraction of gaps that exceeds max_gap_fraction, it is discarded. Set this value to 1 to keep all the sequences.\n\nThe input file can be plaintext (ASCII) or gzip-compressed plaintext (with the extension \".gz\")\n\nSee also remove_duplicate_sequences.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"remove_duplicate_sequences","category":"page"},{"location":"#DCAUtils.remove_duplicate_sequences","page":"Home","title":"DCAUtils.remove_duplicate_sequences","text":"remove_duplicate_sequences(Z::Matrix{Int8}; verbose::Bool = true) -> (Matrix{Int8}, Vector{Int})\n\nTakes a matrix representing a mutiple sequence alignment (see read_fasta_alignment) and returns a new matrix with all duplicated sequences removed. It also returns a vector of column indices with the positions of the unique sequences in the input matrix.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"compute_theta","category":"page"},{"location":"#DCAUtils.compute_theta","page":"Home","title":"DCAUtils.compute_theta","text":"compute_theta(Z::Matrix{Int8}) -> Float64\n\nThis function computes the threshold as used by the :auto setting of compute_weights (but it is more efficient computationally to use the :auto option than to invoke this function and passing the result to compute_weights).\n\nIt computes the mean value ϕ of the similarity fraction between all possible pairs of sequences in Z (the similarity fraction is the number of equal entries divided by the length of the sequences).\n\nThe result is then computed as θ = min(05 01216  ϕ).\n\nThis function can use multiple threads if available.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"compute_weights","category":"page"},{"location":"#DCAUtils.compute_weights","page":"Home","title":"DCAUtils.compute_weights","text":"compute_weights(Z::Matrix{Int8}, [q::Integer,] θ; verbose::Bool = true) -> (Vector{Float64}, Float64)\n\nThis function computes the reweighting vector. It retuns the vector and its sum Meff (the latter represents the number of \"effective\" sequences).\n\nZ is an N  M multiple sequence aLignment (see read_fasta_alignment). N is the length of each sequence and M the number of sequences.\n\nq is the maximum value in the alphabet, if omitted it's computed from maximum(Z).\n\nθ is the distance threshold: for any sequence, the number n of sequences (including itself) that are at normalized distance smaller than θN is counted, and the weight of that sequence is then 1n.\n\nθ can be a real value between 0 and 1, or the symbol :auto, in which case the compute_theta function is used.\n\nThis function can use multiple threads if available.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"compute_dists","category":"page"},{"location":"#DCAUtils.compute_dists","page":"Home","title":"DCAUtils.compute_dists","text":"compute_dists(Z::Matrix{Int8}) -> Matrix{Float64}\n\nThis function computes the matrix of normalized Hamming distances between sequences of the multiple sequence alignment Z (see read_fasta_alignment).\n\nThis function can use multiple threads if available.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"compute_weighted_frequencies","category":"page"},{"location":"#DCAUtils.compute_weighted_frequencies","page":"Home","title":"DCAUtils.compute_weighted_frequencies","text":"compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}) -> (Vector{Float64}, Matrix{Float64})\n\nGiven a multiple sequence alignment matrix Z (see read_fasta_alignment), and a vector of weights (see compute_weights), returns the empirical one- and two-point frequencies P_i and P_ij.\n\nIf Z has size N  M (i.e. M sequences of length N), and its maximum value (the size of the alphabet) is q, the resulting vector P_i has length N (q-1) and contains N blocks (one for each residue position), each block containing the frequencies of the amino-acids, weighted according to W.  The frequency of the last symbol, which usually represents the gap, is omitted and can be recovered by normalization. The resulting matrix P_ij has size N (q-1)  N (q-1) and it also has a block structure, with N  N blocks, one for each pair of residues (the last row and column of each block are omitted and can be recovered by normalization).\n\ncompute_weighted_frequencies(Z::Matrix{Int8}, [q,] θ) -> (Vector{Float64}, Matrix{float64}, Float64, Vector{Float64})\n\nThis form of the function just calls compute_weights with the given values of θ and q and then uses the result to call the version desrcibed above.\n\nBesides returning the one- and two-point frequencies, it also returns the result of compute_weights: the Meff and the reweighting vector.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"add_pseudocount","category":"page"},{"location":"#DCAUtils.add_pseudocount","page":"Home","title":"DCAUtils.add_pseudocount","text":"add_pseudocount(Pi::Vector{Float64}, Pij::Matrix{Float64}, pc::Float64, q::Integer = 21) -> (Vector{Float64}, Matrix{Float64})\n\nThis function takes one- and two-points frequencies (see compute_weighted_frequencies) and returns the corresponding frequencies with a pseudocount pc added.\n\nThe resulting frequencies are the same that would be obtained by mixing the original data with weight (1-pc) with a uniform distribution with weight pc. So pc must be between 0 (returns a copy of the original data) and 1 (returns the frequencies for the uniform distribution).\n\nThe integer q is used to specify the block size of the matrices (each block has size (q-1)(q-1)).\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"compute_DI_gauss","category":"page"},{"location":"#DCAUtils.compute_DI_gauss","page":"Home","title":"DCAUtils.compute_DI_gauss","text":"compute_DI_gauss(J::Matrix{Float64}, C::Matrix{Float64}, q::Integer = 21) -> Matrix{Float64}\n\nCompute the Direct Information matrix, assuming a Gaussian model.\n\nC is the covariance matrix C_ij = P_ij - P_i P_j, and J its inverse.\n\nThe integer q is used to specify the block size of the matrices (each block has size (q-1)(q-1)). The result has one entry per block.\n\nThis function can use multiple threads if available.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"compute_FN","category":"page"},{"location":"#DCAUtils.compute_FN","page":"Home","title":"DCAUtils.compute_FN","text":"compute_FN(J::Matrix{Float64}, q::Integer = 21) -> Matrix{Float64}\n\nCompute the matrix of Frobenius norms.\n\nJ is the inverse of the covariance matrix C_ij = P_ij - P_i P_j.\n\nThe integer q is used to specify the block size of the matrices (each block has size (q-1)(q-1)). The result has one entry per block.\n\nThis function can use multiple threads if available.\n\n\n\n\n\n","category":"function"}]
}