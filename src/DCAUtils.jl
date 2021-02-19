module DCAUtils

export read_fasta_alignment,
       remove_duplicate_sequences,
       compute_theta,
       compute_weights,
       compute_dists,
       compute_weighted_frequencies,
       add_pseudocount,
       compute_DI_gauss,
       compute_FN

include("read_fasta_alignment.jl")
using .ReadFastaAlignment

include("compress_sequencies.jl")
using .CompressSequencies

include("parallel_triu.jl")
using .ParallelTriU

macro hash(x)
    tst = get(ENV, "DCAUTILS_TESTING", "false") == "true"
    return Expr(:call, tst ? :sum : :hash, esc(x))
end

"""
    remove_duplicate_sequences(Z::Matrix{Int8}; verbose::Bool = true) -> (Matrix{Int8}, Vector{Int})

Takes a matrix representing a mutiple sequence alignment (see [`read_fasta_alignment`](@ref))
and returns a new matrix with all duplicated sequences removed. It also returns a vector of column
indices with the positions of the unique sequences in the input matrix.
"""
function remove_duplicate_sequences(Z::Matrix{Int8}; verbose::Bool = true)
    N, M = size(Z)
    hZ = Array{UInt}(undef, M)
    @inbounds for i = 1:M
        hZ[i] = @hash(@view Z[:,i])
    end
    verbose && print("removing duplicate sequences... ")

    ref_seq_ind = Array{Int}(undef, M)
    ref_seq = Dict{UInt,Int}()
    @inbounds for i = 1:M
        ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
    end
    uniqueseqs = collect(values(ref_seq))

    # Check for collisions
    old_collided = trues(M)
    collided = falses(M)
    while true
        fill!(collided, false)
        @inbounds for i = 1:M
            k = ref_seq_ind[i]
            (!old_collided[i] || k == i) && continue
            collided[i] = @views Z[:,i] ≠ Z[:,k]
        end
        any(collided) || break

        # Collect index of first row for each collided hash
        empty!(ref_seq)
        @inbounds for i = 1:M
            collided[i] || continue
            ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
        end
        for v in values(ref_seq)
            push!(uniqueseqs, v)
        end
        old_collided, collided = collided, old_collided
    end
    sort!(uniqueseqs)

    newM = length(uniqueseqs)
    newZ = Z[:,uniqueseqs]
    verbose && println("done: $M -> $newM")
    return newZ, uniqueseqs
end

"""
    compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}) -> (Vector{Float64}, Matrix{Float64})

Given a multiple sequence alignment matrix `Z` (see [`read_fasta_alignment`](@ref)), and a vector
of weights (see [`compute_weights`](@ref)), returns the empirical one- and two-point frequencies
\$P_i\$ and \$P_{ij}\$.

If `Z` has size \$N × M\$ (i.e. \$M\$ sequences of length \$N\$), and its maximum value (the size
of the alphabet) is \$q\$, the resulting vector \$P_i\$ has length \$N (q-1)\$ and contains \$N\$
blocks (one for each residue position), each block containing the frequencies of the amino-acids,
weighted according to `W`.  The frequency of the last symbol, which usually represents the gap, is
omitted and can be recovered by normalization. The resulting matrix \$P_{ij}\$ has size \$N (q-1) ×
N (q-1)\$ and it also has a block structure, with \$N × N\$ blocks, one for each pair of residues
(the last row and column of each block are omitted and can be recovered by normalization).


    compute_weighted_frequencies(Z::Matrix{Int8}, [q,] θ) -> (Vector{Float64}, Matrix{float64}, Float64, Vector{Float64})

This form of the function just calls [`compute_weights`](@ref) with the given values of `θ` and `q`
and then uses the result to call the version desrcibed above.

Besides returning the one- and two-point frequencies, it also returns the result of
`compute_weights`: the `Meff` and the reweighting vector.
"""
function compute_weighted_frequencies end

compute_weighted_frequencies(Z::Matrix{Int8}, θ::Union{Real,Symbol}) = compute_weighted_frequencies(Z, maximum(Z), θ)

function compute_weighted_frequencies(Z::Matrix{Int8}, q::Integer, θ::Union{Real,Symbol})
    W, Meff = compute_weights(Z, q, θ)
    Pi_true, Pij_true = compute_weighted_frequencies(Z, W, Meff)
    return Pi_true, Pij_true, Meff, W
end

compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}) = compute_weighted_frequencies(Z, W, sum(W))

function compute_weighted_frequencies(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64)
    N, M = size(Z)
    q = maximum(Z)
    s = q - 1

    Ns = N * s

    Pij = zeros(Ns, Ns)
    Pi = zeros(Ns)

    ZZ = Vector{Int8}[Z[i,:] for i = 1:N]

    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        for k = 1:M
            a = Zi[k]
            a == q && continue
            Pi[i0 + a] += W[k]
        end
        i0 += s
    end
    Pi ./= Meff

    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        j0 = i0
        for j = i:N
            Zj = ZZ[j]
            for k = 1:M
                a = Zi[k]
                b = Zj[k]
                (a == q || b == q) && continue
                Pij[i0+a, j0+b] += W[k]
            end
            j0 += s
        end
        i0 += s
    end
    for i = 1:Ns
        Pij[i,i] /= Meff
        for j = i+1:Ns
            Pij[i,j] /= Meff
            Pij[j,i] = Pij[i,j]
        end
    end

    return Pi, Pij
end

"""
    add_pseudocount(Pi::Vector{Float64}, Pij::Matrix{Float64}, pc::Float64, q::Integer = 21) -> (Vector{Float64}, Matrix{Float64})

This function takes one- and two-points frequencies (see [`compute_weighted_frequencies`](@ref))
and returns the corresponding frequencies with a pseudocount `pc` added.

The resulting frequencies are the same that would be obtained by mixing the original data with
weight `(1-pc)` with a uniform distribution with weight `pc`. So `pc` must be between 0 (returns a
copy of the original data) and 1 (returns the frequencies for the uniform distribution).

The integer `q` is used to specify the block size of the matrices (each block has size
\$(q-1)×(q-1)\$).
"""
function add_pseudocount(Pi_true::Vector{Float64}, Pij_true::Matrix{Float64}, pc::Float64, q::Integer = 21)
    Nq = length(Pi_true)
    s = q - 1
    Nq % s == 0 || throw(ArgumentError("incompatible length of Pi_true with q"))
    size(Pij_true) == (Nq, Nq) || throw(ArgumentError("incompatible sizes of Pi_true and Pij_true"))
    N = Nq ÷ s

    pcq = pc / q

    Pij = @. (1 - pc) * Pij_true + pcq / q
    Pi = @. (1 - pc) * Pi_true + pcq

    i0 = 0
    for i = 1:N
        xr = i0 .+ (1:s)
        Pij[xr, xr] .= @. (1 - pc) * @view Pij_true[xr, xr]
        for α = 1:s
            x = i0 + α
            Pij[x, x] += pcq
        end
        i0 += s
    end

    return Pi, Pij
end

include("compute_theta.jl")
include("compute_weights.jl")
include("compute_dists.jl")
include("compute_DI_gauss.jl")
include("compute_FN.jl")

end # module
