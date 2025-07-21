

using StaticArrays
using LinearAlgebra
using BenchmarkTools


struct Mesh{TT,NVERTS,NTRIS}
    vertices :: NTuple{NVERTS,SVector{3,TT}}
    connections :: NTuple{NTRIS,NTuple{3,Integer}}

    Mesh{TT}(vertices,connections) where TT = new{TT,length(vertices),length(connections)}(vertices,connections)
end
Base.show(::Mesh{TT,NVERTS,NTRIS}) where {TT,NVERTS,NTRIS} = println("Mesh with $(NVERTS) vertices and $(NTRIS) triangles.")




function toroidal_surface_point(θ,ζ,R₀=5.0,a=1.0)
    x = (R₀ + a * cos(θ)) * cos(ζ)
    y = -(R₀ + a * cos(θ)) * sin(ζ)
    z = a * sin(θ)
    return x, y, z
end





function triangulate_toroidal_surface(θₙ,ζₙ)
    Δθ = π/Float64(θₙ)
    Δζ = π/Float64(ζₙ)

    θ = range(0.0+Δθ,2π+Δθ,θₙ+1)[1:end-1]
    ζ = range(0.0+Δζ,2π+Δζ,ζₙ+1)[1:end-1]

    vertices = Tuple(SVector(toroidal_surface_point(t,z)) for t in θ for z in ζ)

    # Upper and lower triangle shape are
    # (i, j+1) --- (i+1, j+1)
    #   |   \       |
    #   |    \      |
    # (i, j) --- (i+1, j)
    # In linear indexing
    #
    # (i + ζₙ) -- (i + ζₙ + 1)
    #   |   \       |
    # (i)    --   (i + 1)
    #
    #
    # Lower triangle
    # (i, i + ζₙ, i + 1) for i in 1:θₙ
    # Upper triangle
    # (i + 1, i + ζₙ, i + ζₙ + 1)

    ζθ_indices = LinearIndices((ζₙ,θₙ))
    
    lower_tris = ((ζθ_indices[i], ζθ_indices[mod1(i+ζₙ,ζₙ*θₙ)], ζθ_indices[mod1(i+1,ζₙ*θₙ)]) for i in 1:ζₙ*θₙ)
    upper_tris = ((ζθ_indices[mod1(i+1,ζₙ*θₙ)], ζθ_indices[mod1(i+ζₙ,ζₙ*θₙ)], ζθ_indices[mod1(i+1+ζₙ,ζₙ*θₙ)]) for i in 1:ζₙ*θₙ)

    connections =  Iterators.flatten((lower_tris,upper_tris)) |> collect

    return vertices, Tuple(connections)
end





θₙ = 100
ζₙ = 500

verts,conns = triangulate_toroidal_surface(θₙ,ζₙ);

m = Mesh{Float64}(verts,conns);







"""
    (x₁,x₂,x₃) × (y₁,y₂,y₃)
    x₂y₃ - x₃y₂,
    - x₁y₃ + y₁x₃,
    x₁y₂ - x₂y₁
"""
cross(x,y) = (
    x[2]*y[3] - x[3]-y[2],
    y[1]*x[3] - y[3]*x[1],
    x[1]*y[2] - x[2]*y[1]
)


vertex_sign(sign_x,sign_y) = iszero(sign_x) ? sign_x : sign_y
vertex_sign(sign_x,sign_y,sign_z) = vertex_sign(vertex_sign(sign_x,sign_y),sign_z)




"""
Compute the winding number on a single triangle
"""
function winding_number(vertex1,vertex2,vertex3,point)

    # vertex1 - point
    # vertex2 - point
    # vertex3 - point

    # Check the sign of the 
    #   sign(v - pt)
    # if positive 
    signed_vertex1 = vertex_sign((vertex1 - point)...)
    signed_vertex2 = vertex_sign((vertex2 - point)...)
    signed_vertex3 = vertex_sign((vertex3 - point)...)


    check_sign = vertex_sign(signed_vertex1, signed_vertex2, signed_vertex3)
    
    
    sign.(cross(vertex1,point))


    

end






function winding_number(m::Mesh{TT,NVERTS,NTRIS}) where {TT,NVERTS,NTRIS} end



function winding_number(m::Mesh{TT},point::NTuple{3,TT}) where TT
    rx, ry, rz = m.vertices - point
end













"""
Compute the solid angle of a single triangle

Ω/2 = atan.(v₁(v₂×v₃), |v₁||v₂||v₃| + (v₁⋅v₂)|v₃| + (v₂⋅v₃)|v₁|)
"""
function solid_angle(v1,v2,v3,pt)
    

    v1p = v1 - pt
    v2p = v2 - pt
    v3p = v3 - pt

    v1n, v2n, v3n = norm(v1p), norm(v2p), norm(v3p)

    (v1p-pt)*cross(v2p-pt,v3p-pt)

    norm(v1p) * norm(v2p) * norm(v3p) + (v1p*v2p)*v3n + (v1p*v3p)*v2n + (v2p*v3p)*v1n

    return atan(numer,denom)/2π    

end


q = [0.0,0.0,0.0]

solid_angle(m.vertices[m.connections[1][1]],m.vertices[m.connections[1][2]],m.vertices[m.connections[1][3]],q)

