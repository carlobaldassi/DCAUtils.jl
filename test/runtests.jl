module DCAUtilsTests

using DCAUtils
using Test
using JLD2
using LinearAlgebra
using LazyArtifacts

datadir = artifact"test_data"

fastafile(name) = joinpath(datadir, "$name.fasta.gz")

@info "Loading results data"
@load joinpath(datadir, "results.jld2") results

@testset "read FASTA alignment" begin
    for name in ["small", "large"]
        resdict = results[name]
        fn = fastafile(name)
        for (mgf, rd) in keys(resdict["Z"])
            if !rd
                Z = read_fasta_alignment(fn, mgf)
                @test Z == resdict["Z"][(mgf, rd)]
            else
                Z = read_fasta_alignment(fn, mgf)
                Z, _ = remove_duplicate_sequences(Z; verbose=false)
                @test Z == resdict["Z"][(mgf, rd)]
            end
        end
    end
end

@testset "compute θ" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd) in keys(resdict["θ"])
            (mgf, rd) in keys(resdict["Z"]) || continue
            Z = resdict["Z"][(mgf,rd)]
            θ = compute_theta(Z)
            @test θ ≈ resdict["θ"][(mgf,rd)]
            # test slow fallback too
            θ2 = compute_theta([Z[:,i] for i in 1:size(Z,2)], size(Z)...)
            @test θ2 ≈ resdict["θ"][(mgf,rd)]
        end
    end
end

@testset "compute dists" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd) in keys(resdict["D"])
            (mgf, rd) in keys(resdict["Z"]) || continue
            Z = resdict["Z"][(mgf,rd)]
            D = compute_dists(Z)
            @test D ≈ resdict["D"][(mgf,rd)]
            # test slow fallback too
            D2 = compute_dists([Z[:,i] for i in 1:size(Z,2)], size(Z)...)
            @test D2 ≈ resdict["D"][(mgf,rd)]
        end
    end
end

@testset "compute weights" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd, θ) in keys(resdict["W"])
            (mgf, rd) in keys(resdict["Z"]) || continue
            Z = resdict["Z"][(mgf,rd)]
            W, Meff = compute_weights(Z, θ; verbose = false)
            @test Meff ≈ sum(W)
            @test W ≈ resdict["W"][(mgf, rd, θ)]
        end
    end
end

@testset "compute weighted frequencies" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd, θ) in keys(resdict["Pi_true"])
            (mgf, rd) in keys(resdict["Z"]) || continue
            (mgf, rd, θ) in keys(resdict["W"]) || continue
            (mgf, rd, θ) in keys(resdict["Pij_true"]) || continue

            Z = resdict["Z"][(mgf,rd)]
            W = resdict["W"][(mgf, rd, θ)]
            Pi_true, Pij_true = compute_weighted_frequencies(Z, W)
            @test Pi_true ≈ resdict["Pi_true"][(mgf, rd, θ)]
            @test Pij_true ≈ resdict["Pij_true"][(mgf, rd, θ)]
        end
    end
end

@testset "add pseudocount" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd, θ, pc) in keys(resdict["Pi_pc"])
            (mgf, rd, θ) in keys(resdict["Pi_true"]) || continue
            (mgf, rd, θ) in keys(resdict["Pij_true"]) || continue
            (mgf, rd, θ, pc) in keys(resdict["Pij_pc"]) || continue
            Pi_true = resdict["Pi_true"][(mgf, rd, θ)]
            Pij_true = resdict["Pij_true"][(mgf, rd, θ)]
            Pi_pc, Pij_pc = add_pseudocount(Pi_true, Pij_true, pc)

            @test Pi_pc ≈ resdict["Pi_pc"][(mgf, rd, θ, pc)]
            @test Pij_pc ≈ resdict["Pij_pc"][(mgf, rd, θ, pc)]
        end
    end
end

@testset "compute direct information" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd, θ, pc) in keys(resdict["DI"])
            (mgf, rd, θ, pc) in keys(resdict["Pi_pc"]) || continue
            (mgf, rd, θ, pc) in keys(resdict["Pij_pc"]) || continue
            Pi_pc = resdict["Pi_pc"][(mgf, rd, θ, pc)]
            Pij_pc = resdict["Pij_pc"][(mgf, rd, θ, pc)]
            C = Pij_pc - Pi_pc * Pi_pc'
            mJ = inv(cholesky(C))
            DI = compute_DI(mJ, C)
            @test DI ≈ resdict["DI"][(mgf, rd, θ, pc)]
        end
    end
end

@testset "compute Frobenius norm" begin
    for name in ["small", "large"]
        resdict = results[name]
        for (mgf, rd, θ, pc) in keys(resdict["FN"])
            (mgf, rd, θ, pc) in keys(resdict["Pi_pc"]) || continue
            (mgf, rd, θ, pc) in keys(resdict["Pij_pc"]) || continue
            Pi_pc = resdict["Pi_pc"][(mgf, rd, θ, pc)]
            Pij_pc = resdict["Pij_pc"][(mgf, rd, θ, pc)]
            C = Pij_pc - Pi_pc * Pi_pc'
            mJ = inv(cholesky(C))
            FN = compute_FN(mJ)
            @test FN ≈ resdict["FN"][(mgf, rd, θ, pc)]
        end
    end
end

end # module
