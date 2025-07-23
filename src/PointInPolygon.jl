module PointInPolygon

using StaticArrays
using LinearAlgebra

export Mesh
export solid_angle, solid_angle!
export winding_number, winding_number!


struct Mesh{TT,NVERTS,NTRIS}
    vertices::NTuple{NVERTS,SVector{3,TT}}
    connections::NTuple{NTRIS,SVector{3,Int}}

    Mesh{TT}(vertices, connections) where {TT} = new{TT,length(vertices),length(connections)}(vertices, connections)
end

Base.show(::Mesh{TT,NVERTS,NTRIS}) where {TT,NVERTS,NTRIS} = println("Mesh with $(NVERTS) vertices and $(NTRIS) triangles.")


"""
    (x₁,x₂,x₃) × (y₁,y₂,y₃)
    x₂y₃ - x₃y₂,
    - x₁y₃ + y₁x₃,
    x₁y₂ - x₂y₁
"""
cross(x, y) = SVector(
    x[2] * y[3] - x[3] * y[2],
    y[1] * x[3] - y[3] * x[1],
    x[1] * y[2] - x[2] * y[1]
)
"""
    (x-z)×(y-z)
"""
crossm(x, y, z) = SVector(x - z, y - z)


"""
Compute the solid angle of a single triangle

Ω/2 = atan(v₁(v₂×v₃), |v₁||v₂||v₃| + (v₁⋅v₂)|v₃| + (v₂⋅v₃)|v₁|)
"""
function solid_angle(v1, v2, v3, pt)

    v1p = v1 - pt
    v2p = v2 - pt
    v3p = v3 - pt

    v1n, v2n, v3n = norm(v1p), norm(v2p), norm(v3p)

    numer = dot(v1p, cross(v2p, v3p))

    denom = v1n * v2n * v3n + dot(v1p, v2p) * v3n + dot(v1p, v3p) * v2n + dot(v2p, v3p) * v1n

    return atan(numer, denom) / (2π)
end
function solid_angle(mesh::Mesh{TT,NVERTS,NTRIS}, point::AbstractArray{TT}) where {TT,NVERTS,NTRIS}
    Ω = mapreduce(connection -> solid_angle(mesh.vertices[connection[1]],
            mesh.vertices[connection[2]],
            mesh.vertices[connection[3]],
            point), +, mesh.connections)
    return Ω
end
solid_angle(mesh::Mesh, points::NTuple) = map(pt -> solid_angle(mesh, pt), points)

function solid_angle!(result, mesh::Mesh, points::NTuple)
    for i in eachindex(result)
        result[i] = solid_angle(mesh, points[i])
    end
end


vertex_sign(x, y) = iszero(x) ? sign(y) : sign(x)
vertex_sign(x, y, z) = vertex_sign(vertex_sign(x, y), z)
vertex_sign(x) = vertex_sign(x...)

edge_sign(v₁, v₂) = (v₁[2] * v₂[1] - v₁[1] * v₂[2], v₁[3] * v₂[1] - v₁[1] * v₂[3], v₁[3] * v₂[2] - v₁[2] * v₂[3])

triangle_area(v₁, v₂, v₃) = (
    (v₁[1] * v₂[2] - v₁[2] * v₂[1]) * v₃[3] +
    (v₂[1] * v₃[2] - v₂[2] * v₃[1]) * v₁[2] +
    (v₃[1] * v₁[2] - v₃[2] * v₂[1]) * v₂[2]
)

"""
Compute the winding number on a single triangle
"""
function winding_number(v1, v2, v3, point::AbstractVector{TT}) where {TT}

    v1p = v1 - point
    v2p = v2 - point
    v3p = v3 - point

    v1_sign = vertex_sign(v1p)
    v2_sign = vertex_sign(v2p)
    v3_sign = vertex_sign(v3p)

    check_faces = zero(TT)

    if v1_sign != v2_sign
        check_faces += vertex_sign(edge_sign(v1p, v2p))
    end
    if v2_sign != v3_sign
        check_faces += vertex_sign(edge_sign(v2p, v3p))
    end
    if v3_sign != v1_sign
        check_faces += vertex_sign(edge_sign(v3p, v1p))
    end


    winding_number_contribution = zero(TT)
    if !iszero(check_faces)
        winding_number_contribution += sign(triangle_area(v1p, v2p, v3p))
    end
    return winding_number_contribution
end
function winding_number(mesh::Mesh{TT,NVERTS,NTRIS}, point::AbstractArray{TT}) where {TT,NVERTS,NTRIS}
    wₙ = mapreduce(connection -> winding_number(
            mesh.vertices[connection[1]],
            mesh.vertices[connection[2]],
            mesh.vertices[connection[3]],
            point), +, mesh.connections)
    return fld(wₙ, 2)
end
winding_number(mesh::Mesh, points::NTuple) = map(pt -> winding_number(mesh, pt), points)

function winding_number!(result, mesh::Mesh, points::NTuple)
    for i in eachindex(result)
        result[i] = winding_number(mesh, points[i])
    end
end

end
