module Example001_Solvers

# under development

using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearSolve
using ILUZero
using LinearAlgebra



function main(; n = 10, Plotter = nothing, kwargs...)

    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = VoronoiFVM.Grid(X, Y)

    eps = 1.0e-2

    function reaction(f, u, node)
        f[1] = u[1]^2
    end

    function flux(f, u, edge)
        f[1] = eps * (u[1, 1]^2 - u[1, 2]^2)
    end

    function source(f, node)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20.0 * (x1^2 + x2^2))
    end

    function storage(f, u, node)
        f[1] = u[1]
    end

    function bcondition(f, u, node)
        boundary_dirichlet!(
            f,
            u,
            node,
            species = 1,
            region = 2,
            value = ramp(node.time, dt = (0, 0.1), du = (0, 1)),
        )
        boundary_dirichlet!(
            f,
            u,
            node,
            species = 1,
            region = 4,
            value = ramp(node.time, dt = (0, 0.1), du = (0, 1)),
        )
    end

    sys =
        VoronoiFVM.System(grid; reaction, flux, source, storage, bcondition, species = [1])

    @info "KLU:"

    klu_sol = solve(sys; inival = 0.5, method_linear = KLUFactorization(), kwargs...)

    @info "UMFPACK:"
    umf_sol = solve(sys; inival = 0.5, method_linear = UMFPACKFactorization(), kwargs...)

    @info "Sparspak:"
    spk_sol = solve(sys; inival = 0.5, method_linear = SparspakFactorization(), kwargs...)

    @info "Krylov:"

    kry_sol = solve(
        sys;
        inival = 0.5,
        method_linear = KrylovJL_BICGSTAB(),
        precon_linear = ILUZero.ilu0,
        kwargs...,
    )


    @info "Krylov2:"
    kry2_sol = solve(
        sys;
        inival = 0.5,
        method_linear = KrylovJL_BICGSTAB(),
        precon_linear = A -> VoronoiFVM.factorization(A, SparspakFactorization()),
        kwargs...,
    )

    klu_sol ≈ umf_sol && spk_sol ≈ umf_sol
    norm(kry_sol - umf_sol, Inf)
end

function test()
    main()
    true
end
end
