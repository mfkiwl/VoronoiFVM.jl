module test_cachesol
using ExtendableGrids
using VoronoiFVM
using LinearSolve
using ExtendableSparse
using LinearAlgebra
using SparseArrays
using Test

function flux(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
    return nothing
end

function bcondition(f, u, bnode, data)
    boundary_dirichlet!(f, u, bnode; species = 1, region = 1, value = 0.0)
    boundary_dirichlet!(f, u, bnode; species = 1, region = 2, value = 1.0)
    return nothing
end

function setup_linear_solve(n)
    # define system, initialize system matrix
    X = linspace(0.0, 1.0, n)
    system = VoronoiFVM.System(X; flux, bcondition, species = [1])
    state = VoronoiFVM.SystemState(system)
    F = unknowns(system; inival = 0.0)
    Tv = eltype(F)
    # assemble system matrix A s.th. the solution of the system is A*F=RHS
    # for RHS correctly assembled and the state.residual = A*F - RHS;
    # so for F=0, we have in state.residual = -RHS the correctly assembled RHS
    # such that X solves the discretized system A*X = RHS
    VoronoiFVM.eval_and_assemble(
        state.system,
        F, # current solution where A*F is evaluated
        F, # previous time step solution, doesn't do anything here
        state.residual, # A*F - RHS = A*0 - RHS = -RHS
        state.matrix, # matrix A
        state.dudp,
        0.0, # time value
        Inf, # tstep value
        0.0, # embedding parameter
        state.data,
        Tv[] # system parameters
    )
    F .= -state.residual

    ## working version: pull out the sparse matrix and
    #  vectorize RHS, dispatch to LinearSolve
    #flush!(state.matrix)
    #mtx = state.matrix.cscmatrix
    #p = LinearProblem(mtx, VoronoiFVM.dofs(F))
    #linear_cache = init(p, LinearSolve.UMFPACKFactorization())
    #linear_cache.b = VoronoiFVM.dofs(F)

    ## broken version as of LinearSolve v. 3.7.0: initiate problem
    #  and RHS as described in the docs of ExtendableSparse and LinearSolve
    mtx = state.matrix
    p = LinearProblem(mtx, VoronoiFVM.dofs(F))
    linear_cache = init(p, LinearSolve.UMFPACKFactorization())
    linear_cache.b = VoronoiFVM.dofs(F)

    return X, linear_cache, state, mtx
end

function solve_problem(cache)
    sol = LinearSolve.solve!(cache)
    return sol
end

function main(; n = 10)
    X, cache, _ = setup_linear_solve(n)
    approx_sol = solve_problem(cache)
    true_sol = map((x) -> (x), X)
    @test norm(approx_sol .- true_sol, Inf) ≤ 3.0e-16
    return nothing
end

function runtests()
    main()
    return nothing
end
end
