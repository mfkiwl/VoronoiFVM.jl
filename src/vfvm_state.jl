"""
$(TYPEDEF)

Structure holding state information for finite volume system.


Type parameters:
- Tv: element type of solution vectors and matrix
- TMatrix:  matrix type
- TSolArray: type of solution vector: (Matrix or SparseMatrixCSC)
- TData: type of user data

Type fields:
$(TYPEDFIELDS)
"""
mutable struct SystemState{Tv, TMatrix<:AbstractMatrix{Tv}, TSolArray<:AbstractMatrix{Tv}, TData}

    """
    Related finite volume system
    """
    system::VoronoiFVM.System

    """
    Parameter data 
    """
    data::TData

    """
    Solution vector
    """
    solution::TSolArray

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::TMatrix

    """
    Parameter derivative (vector of solution arrays)
    """
    dudp::Vector{TSolArray}
    
    """
    Vector holding Newton update
    """
    update::TSolArray

    """
    Vector holding Newton residual
    """
    residual::TSolArray

    """
    Linear solver cache
    """
    linear_cache::Union{Nothing, LinearSolve.LinearCache}

    """
    Hash value of latest unknowns vector the assembly was called with
    (used by differential equation interface)
    """
    uhash::UInt64

    """
    History record for solution process
    (used by differential equation interface)
    """
    history::Union{DiffEqHistory, Nothing}

end


"""
    SystemState(Tv, system; data=system.physics.data, matrixtype=system.matrixtype)

Create state information for finite volume system.

Arguments:
- `Tv`: value type of unknowns, matrix
- `system`: Finite volume system

Keyword arguments:
- `data`: User data. Default: `data(system)`
- `matrixtype`. Default: `system.matrixtype`
"""
function SystemState(::Type{Tu}, system::AbstractSystem{Tv, Tc, Ti, Tm};
                     data=system.physics.data,
                     matrixtype=system.matrixtype) where {Tu,Tv,Tc, Ti, Tm}
    _complete!(system)

    nspec = size(system.node_dof, 1)
    n = num_dof(system)

    matrixtype = system.matrixtype

    if matrixtype == :default
        if !isdensesystem(system) || dim_grid(system.grid)>1
            matrixtype = :sparse
        else
            if nspec == 1
                matrixtype = :tridiagonal
            else
                matrixtype = :banded
            end
        end
    end

    if matrixtype == :tridiagonal
        matrix = Tridiagonal(zeros(Tu, n - 1), zeros(Tu, n), zeros(Tu, n - 1))
    elseif matrixtype == :banded
        matrix = BandedMatrix{Tu}(Zeros(n, n), (2 * nspec - 1, 2 * nspec - 1))
        # elseif matrixtype==:multidiagonal
        #     system.matrix=mdzeros(Tv,n,n,[-1,0,1]; blocksize=nspec)
    else # :sparse
        if num_partitions(system.grid) == 1
            matrix = ExtendableSparseMatrixCSC{Tu, Tm}(n, n)
        else
            matrix = MTExtendableSparseMatrixCSC{Tu, Tm}(n, n, num_partitions(system.grid))
        end
    end

    solution = unknowns(Tu, system)
    residual = unknowns(Tu, system)
    update = unknowns(Tu, system)
    dudp = [unknowns(Tu, system) for i = 1:(system.num_parameters)]
    SystemState(system, data, solution, matrix, dudp, residual, update, nothing, zero(UInt64),  nothing)
end


"""
    SystemState(system; kwargs...)

Shortcut for creating state with value type defined by `Tv` type parameter of system
"""
SystemState(system::AbstractSystem{Tv,Tc,Ti,Tm}; kwargs...) where {Tv,Tc,Ti,Tm} =SystemState(Tv, system; kwargs...) 


