"""
   $(TYPEDEF)

Abstract type for geometry items (node,bnode,edge, bedge)
"""
abstract type AbstractGeometryItem{Tc<:Number, Ti <:Integer} end


time(item::AbstractGeometryItem)=item.time
embedparam(item::AbstractGeometryItem)=item.embedparam
region(item::AbstractGeometryItem)=item.region


"""
   $(TYPEDEF)

Abstract type for nodes. 

`node[idim]` gives the the corresponding coordinate.
"""
abstract type AbstractNode{Tc<:Number, Ti <:Integer}  <: AbstractGeometryItem{Tc, Ti} end
Base.size(node::AbstractNode)=(size(node.coord)[1],)
Base.getindex(node::AbstractNode, idim)=@inbounds  node.coord[idim,node.index]


"""
    $(TYPEDEF)

Abstract type for data on nodes.
`u[ispec]` accesses value of species at this node.
"""
abstract type AbstractNodeData{Tv<: Number} <: AbstractVector{Tv} end
Base.size(u::AbstractNodeData)=(u.nspec,1)
Base.getindex(u::AbstractNodeData,i)=@inbounds u.val[i]
Base.setindex!(f::AbstractNodeData,v,i)=@inbounds f.val[i]=v


struct DParameters{Tv<:Number} <: AbstractVector{Tv}
    val::Vector{Tv}
    offset::Int32
end

Base.size(p::DParameters)=(length(p.val)-p.offset,1)
Base.getindex(p::DParameters,i)=@inbounds p.val[p.offset+i]

function parameters(u::AbstractNodeData{Tv}) where {Tv<:Number}
    DParameters(u.val,u.nspec)
end



"""
   $(TYPEDEF)

Abstract type for edges 

`edge[idim,inode]` gives coordinate of node.
"""
abstract type AbstractEdge{Tv<:Number, Ti <:Integer}  <: AbstractGeometryItem{Tv, Ti} end
Base.size(edge::AbstractEdge)=(size(edge.coord)[1],2)
Base.getindex(edge::AbstractEdge, idim,inode)=@inbounds  edge.coord[idim,edge.node[inode]]

"""
    $(TYPEDEF)

Abstract type for data on edges.
`u[ispec,inode]` accesses value of species at corresponding node.
"""
abstract type AbstractEdgeData{Tv<: Number} <: AbstractMatrix{Tv} end
Base.size(u::AbstractEdgeData)=(u.n1,2)
Base.getindex(u::AbstractEdgeData,i,j)=@inbounds u.val[(j-1)*u.n1+i]

function parameters(u::AbstractEdgeData{Tv}) where {Tv<:Number}
    DParameters(u.val,2*u.n1)
end


##################################################################
"""
$(TYPEDEF)

Structure holding local node information.

$(TYPEDFIELDS)
"""
mutable struct Node{Tc,Ti} <: AbstractNode{Tc, Ti} 

    """
    Index in grid

    """
    index::Ti
    """
    Inner region number
    """
    region::Ti

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Number of discretization cell the node is invoked from
    """
    icell::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tc}


    """
    Grid cell nodes
    """
    cellnodes::Array{Ti,2}

    
    """
    Grid cell regions
    """
    cellregions::Vector{Ti}

    """
    System time
    """
    time::Float64
    
    """
    Current value of embedding parameter
    """
    embedparam::Float64
    
    Node{Tc,Ti}(sys::AbstractSystem{Tv,Tc,Ti,Tm}) where {Tv,Tc,Ti,Tm}  =new(zero(Ti),0,
                                                                num_species(sys),0,
                                                                coordinates(sys.grid),
                                                                sys.grid[CellNodes],
                                                                sys.grid[CellRegions],
                                                                0,0
                                                                )
end

Node(sys::AbstractSystem{Tv,Tc,Ti,Tm}) where {Tv,Tc,Ti,Tm}=Node{Tc,Ti}(sys)

@inline function _fill!(node::Node,inode,icell)
    node.region=node.cellregions[icell]
    node.index=node.cellnodes[inode,icell]
    node.icell=icell
    nothing
end


"""
    $(TYPEDEF)

Unknown data on node. 
"""
struct NodeUnknowns{Tv,Tc,Ti} <:AbstractNodeData{Tv} 
    val::Vector{Tv}
    nspec::Ti
    geom::Node{Tc,Ti}
end

@inline unknowns(node::Node{Tc,Ti},u::AbstractVector{Tv}) where {Tv,Tc,Ti} = NodeUnknowns{Tv,Tc,Ti}(u,node.nspec,node)

"""
    $(TYPEDEF)

RHS data on node. 
"""
struct NodeRHS{Tv,Tc,Ti} <:AbstractNodeData{Tv}
    val::Vector{Tv}
    nspec::Ti
    geom::Node{Tc,Ti}
end

@inline rhs(node::Node{Tc,Ti}, f::AbstractVector{Tv}) where {Tv,Tc,Ti} = NodeRHS{Tv,Tc,Ti}(f,node.nspec,node)


##################################################################
"""
$(TYPEDEF)

Structure holding local boundary  node information.

$(TYPEDFIELDS)
"""
mutable struct BNode{Tv, Tc, Ti} <: AbstractNode{Tc, Ti}

    """
    Index in grid
    """
    index::Ti
    
    """
    BFace number it is called from
    """
    ibface::Ti

    
    """
    local node number
    """
    ibnode::Ti

    """
    Boundary region number
    """
    region::Ti


    cellregions::Vector{Ti}

    """
    Number of species defined in node
    """
    nspec::Ti

    """
    Grid coordinates
    """
    coord::Matrix{Tc}


    bfacenodes::Array{Ti,2}

    bfaceregions::Vector{Ti}

    allcellregions::Vector{Ti}
    
    bfacecells::ExtendableGrids.Adjacency{Ti}
    
    Dirichlet::Tv

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64

    dirichlet_value::Vector{Tv}
    
    BNode{Tv, Tc,Ti}(sys::AbstractSystem{Tv,Tc,Ti,Tm}) where {Tv,Tc,Ti,Tm}  =new(0,0,0,0,zeros(Ti,2),
                                                                                 num_species(sys),
                                                                                 coordinates(sys.grid),
                                                                                 sys.grid[BFaceNodes],
                                                                                 sys.grid[BFaceRegions],
                                                                                 sys.grid[CellRegions],
                                                                                 sys.grid[BFaceCells],
                                                                                 Dirichlet,0,0,
                                                                                 zeros(Tv,num_species(sys))
                                                                                 )
end
BNode(sys::AbstractSystem{Tv,Tc,Ti,Tm}) where {Tv,Tc,Ti,Tm}=BNode{Tv,Tc,Ti}(sys)


@inline function _fill0!(node::BNode,ibnode,ibface)
    node.ibface=ibface
    node.ibnode=ibnode
    node.region=node.bfaceregions[ibface]
    node.index=node.bfacenodes[ibnode,ibface]
    nothing
end


@inline function _fill!(node::BNode,ibnode,ibface)
    _fill0!(node,ibnode,ibface)
    node.cellregions[1]=0
    node.cellregions[2]=0
    for i=1:num_targets(node.bfacecells,ibface)
        icell=node.bfacecells[i,ibface]
        node.cellregions[i]=node.allcellregions[icell]
    end
end




struct BNodeUnknowns{Tv,Tc,Ti} <:AbstractNodeData{Tv} 
    val::Vector{Tv}
    nspec::Ti
    geom::BNode{Tv,Tc,Ti}
end

@inline unknowns(bnode::BNode{Tv,Tc,Ti},u::AbstractVector{Tv}) where {Tv,Tc,Ti} = BNodeUnknowns{Tv,Tc,Ti}(u,bnode.nspec,bnode)


struct BNodeRHS{Tv,Tc,Ti} <:AbstractNodeData{Tv} 
    val::Vector{Tv}
    nspec::Ti
    geom::BNode{Tv,Tc,Ti}
end

@inline rhs(bnode::BNode{Tv,Tc,Ti}, f::AbstractVector{Tv}) where {Tv,Tc,Ti} = BNodeRHS{Tv,Tc,Ti}(f,bnode.nspec,bnode)




##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct Edge{Tc,Ti}  <: AbstractEdge{Tc, Ti}

    """
    Index in grid
    """
    index::Ti

    """
    Index 
    """
    node::Vector{Ti}

    """
    Inner region number corresponding to edge
    """
    region::Ti

    """
    Number of species defined in edge
    """
    nspec::Ti

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Ti
    
    """
    Grid coordinates
    """
    coord::Matrix{Tc}

    
    cellx::Array{Ti,2}
    edgenodes::Array{Ti,2}
    cellregions::Vector{Ti}
    has_celledges::Bool

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64
    
    Edge{Tc,Ti}(::Nothing) where {Tc,Ti}  =new()
end


function  Edge(sys::AbstractSystem{Tv,Tc,Ti,Tm}) where {Tv,Tc,Ti,Tm} 
    edge=Edge{Tc,Ti}(nothing)

    edge.index=0
    edge.node=[0,0]
    edge.region=0
    edge.nspec=num_species(sys)
    edge.icell=0
    edge.coord=coordinates(sys.grid)
    geom=sys.grid[CellGeometries][1]
    if haskey(sys.grid,CellEdges)
        edge.cellx=sys.grid[CellEdges]
        edge.edgenodes=sys.grid[EdgeNodes]
        edge.has_celledges=true
    else
        edge.cellx=sys.grid[CellNodes]
        edge.edgenodes=local_celledgenodes(geom)
        edge.has_celledges=false
    end
    edge.cellregions=sys.grid[CellRegions]
    edge.time=0
    edge.embedparam=0
    edge
end


@inline function _fill!(edge::Edge,iedge,icell)
    if edge.has_celledges #  cellx==celledges, edgenodes==global_edgenodes
        # If we work with projections of fluxes onto edges,
        # we need to ensure that the edges are accessed with the
        # same orientation without regard of the orientation induced
        # by local cell numbering
        edge.index=edge.cellx[iedge,icell]
        edge.node[1]=edge.edgenodes[1,edge.index]
        edge.node[2]=edge.edgenodes[2,edge.index]
    else # cx==cellnodes, edgenodes== local_edgenodes
        edge.index=0
        edge.node[1]=edge.cellx[edge.edgenodes[1,iedge],icell]
        edge.node[2]=edge.cellx[edge.edgenodes[2,iedge],icell]
    end
    edge.region=edge.cellregions[icell]
    edge.icell=icell
    nothing
end


struct EdgeUnknowns{Tv,Tc,Ti} <:AbstractEdgeData{Tv} 
    val::Vector{Tv}
    n1::Ti
    geom::Edge{Tc,Ti}
end

@inline unknowns(edge::Edge{Tc,Ti},u::AbstractVector{Tv}) where {Tv,Tc,Ti} = EdgeUnknowns{Tv,Tc,Ti}(u,edge.nspec,edge)


struct EdgeRHS{Tv,Tc,Ti} <:AbstractNodeData{Tv} 
    val::Vector{Tv}
    nspec::Ti
    geom::Edge{Tc,Ti}
end

@inline rhs(edge::Edge{Tc,Ti}, f::AbstractVector{Tv}) where {Tv,Tc,Ti} = EdgeRHS{Tv,Tc,Ti}(f,edge.nspec,edge)








##################################################################
"""
$(TYPEDEF)

Structure holding local edge information.

$(TYPEDFIELDS)
"""
mutable struct BEdge{Tc,Ti}  <: AbstractEdge{Tc, Ti}

    """
    Index in grid
    """
    index::Ti

    """
    Index 
    """
    node::Vector{Ti}

    """
    Inner region number corresponding to edge
    """
    region::Ti

    """
    Number of species defined in edge
    """
    nspec::Ti

    """
    Number of discretization cell the edge is invoked from
    """
    icell::Ti
    
    """
    Grid coordinates
    """
    coord::Matrix{Tc}

    bedgenodes::Array{Ti,2}
    bfaceedges::Array{Ti,2}
    bfaceregions::Vector{Ti}

    """
    System time
    """
    time::Float64

    """
    Current value of embedding parameter
    """
    embedparam::Float64

    
    BEdge{Tc,Ti}(::Nothing) where {Tc,Ti}  =new()
end


function  BEdge(sys::AbstractSystem{Tv,Tc,Ti,Tm}) where {Tv,Tc,Ti,Tm} 
    bedge=BEdge{Tc,Ti}(nothing)

    bedge.index=0
    bedge.node=[0,0]
    bedge.region=0
    bedge.nspec=num_species(sys)
    bedge.icell=0
    bedge.coord=coordinates(sys.grid)
       
    bedge.bedgenodes=sys.grid[BEdgeNodes]
    bedge.bfaceedges=sys.grid[BFaceEdges]
    bedge.bfaceregions=sys.grid[BFaceRegions]
    bedge.time=0
    bedge.embedparam=0
    bedge
end


@inline function _fill!(bedge::BEdge,ibedge,ibface)

    bedge.index   = bedge.bfaceedges[ibedge, ibface]
    bedge.node[1] = bedge.bedgenodes[1, bedge.index]
    bedge.node[2] = bedge.bedgenodes[2, bedge.index]
    
    bedge.region = bedge.bfaceregions[ibface]
    bedge.icell  = ibface

    nothing
end

struct BEdgeUnknowns{Tv,Tc,Ti} <:AbstractEdgeData{Tv} 
    val::Vector{Tv}
    n1::Ti
    geom::BEdge{Tc,Ti}
end

@inline unknowns(edge::BEdge{Tc,Ti},u::AbstractVector{Tv}) where {Tv,Tc,Ti} = BEdgeUnknowns{Tv,Tc,Ti}(u,edge.nspec,edge)


struct BEdgeRHS{Tv,Tc,Ti} <:AbstractNodeData{Tv} 
    val::Vector{Tv}
    nspec::Ti
    geom::BEdge{Tc,Ti}
end

@inline rhs(edge::BEdge{Tc,Ti}, f::AbstractVector{Tv}) where {Tv,Tc,Ti}= BEdgeRHS{Tv,Tc,Ti}(f,edge.nspec,edge)


##################################################################
"""
$(TYPEDSIGNATURES)

Return number of species for edge
"""
@inline num_species(edge::AbstractEdge)=edge.nspec


##################################################################
"""
$(TYPEDSIGNATURES)
   
Calculate the length of an edge. 
"""
function meas(edge::AbstractEdge)
    l=0.0
    for i=1:size(edge.coord)[1]
        d=edge.coord[i,edge.node[1]]-edge.coord[i,edge.node[2]]
        l=l+d*d
    end
    return sqrt(l)
end

edgelength(edge::AbstractEdge)=meas(edge)


function project(edge::Edge,vec)
    vh=0.0
    for i=1:size(edge.coord)[1]
        vh+=(edge.coord[i,edge.node[2]]-edge.coord[i,edge.node[1]])*vec[i]
    end
    return vh
end






###############################################################
# Deprecation warnings here ?
"""
$(TYPEDEF)

Wrapper struct for viewing unknowns passed to callback functions
    
$(TYPEDFIELDS)
"""
struct VectorUnknowns{Tv} <:AbstractVector{Tv} 
    val::Vector{Tv}
    n::Int64
    offset::Int64
end


"""
$(TYPEDSIGNATURES)

Construct vector unknowns on edge.
"""
unknowns(edge::AbstractEdge, u::AbstractVector{Tv},i) where Tv = VectorUnknowns{Tv}(u,edge.nspec,(i-1)*(edge.nspec))
Base.getindex(u::VectorUnknowns,i)=@inbounds u.val[u.offset+i]
Base.size(u::VectorUnknowns)=(u.n,)



# For backward compatibility
unknowns(edge,u::AbstractEdgeData)=u
# For backward compatibility
unknowns(edge::Edge, u::AbstractEdgeData{Tv},i) where Tv = VectorUnknowns{Tv}(u.val,edge.nspec,(i-1)*(edge.nspec))

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
viewK(edge::AbstractEdge,u)=unknowns(edge,u,1)


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
viewL(edge::AbstractEdge,u)=unknowns(edge,u,2)

