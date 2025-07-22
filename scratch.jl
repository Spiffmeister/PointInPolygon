

using StaticArrays
using LinearAlgebra
using BenchmarkTools


struct Mesh{TT,NVERTS,NTRIS}
    vertices :: NTuple{NVERTS,SVector{3,TT}}
    connections :: NTuple{NTRIS,NTuple{3,Int}}

    Mesh{TT}(vertices,connections) where TT = new{TT,length(vertices),length(connections)}(vertices,connections)
end
Base.show(io::IO,m::Mesh{TT,NVERTS,NTRIS}) where {TT,NVERTS,NTRIS} = println("Mesh with $(NVERTS) vertices and $(NTRIS) triangles.")




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

    lower_triangle = ((ζθ_indices[i,j], ζθ_indices[i,mod1(j+1,θₙ)], ζθ_indices[mod1(i+1,ζₙ),j]) for i in 1:ζₙ for j in 1:θₙ)
    upper_triangle = ((ζθ_indices[mod1(i+1,ζₙ),j], ζθ_indices[i,mod1(j+1,θₙ)], ζθ_indices[mod1(i+1,ζₙ),mod1(j+1,θₙ)]) for i in 1:ζₙ for j in 1:θₙ)
    connections = Iterators.flatten(zip(lower_triangle,upper_triangle)) |> collect

    return vertices, Tuple(connections)
end










"""
    (x₁,x₂,x₃) × (y₁,y₂,y₃)
    x₂y₃ - x₃y₂,
    - x₁y₃ + y₁x₃,
    x₁y₂ - x₂y₁
"""
cross(x,y) = (
    x[2]*y[3] - x[3]*y[2],
    y[1]*x[3] - y[3]*x[1],
    x[1]*y[2] - x[2]*y[1]
)
"""
    (x-z)×(y-z)
"""
crossm(x,y,z) = (x-z, y-z)


vertex_sign(sign_x,sign_y) = iszero(sign_x) ? sign_y : sign_x
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

Ω/2 = atan(v₁(v₂×v₃), |v₁||v₂||v₃| + (v₁⋅v₂)|v₃| + (v₂⋅v₃)|v₁|)
"""
function solid_angle(v1,v2,v3,pt)
    
    v1p = v1 - pt
    v2p = v2 - pt
    v3p = v3 - pt

    v1n, v2n, v3n = norm(v1p), norm(v2p), norm(v3p)

    numer = dot(v1p,cross(v2p,v3p))
    
    denom = v1n*v2n*v3n + dot(v1p,v2p)*v3n + dot(v1p,v3p)*v2n + dot(v2p,v3p)*v1n
    
    return atan(numer,denom)/(2π)
end
function solid_angle(mesh::Mesh{TT,NVERTS,NTRIS},pt::AbstractArray{TT}) where {TT,NVERTS,NTRIS}
    Ω = TT(0)
    for connection in mesh.connections
        Ω += solid_angle(mesh.vertices[connection[1]],
            mesh.vertices[connection[2]],
            mesh.vertices[connection[3]],
            pt)
    end
    return Ω
end
solid_angle_test(mesh::Mesh,points::NTuple) = map(pt -> solid_angle(mesh,pt), points)







θₙ = 100
ζₙ = 500

verts,conns = triangulate_toroidal_surface(θₙ,ζₙ);

m = Mesh{Float64}(verts,conns);

d = rand(3);


@benchmark solid_angle($m,$d)






"""
    Florians test particles
"""

testpt() = toroidal_surface_point(rand()*2π, rand()*2π, 5.1, rand()*0.27 + 0.82)
florians_test_points(ntest=400) = Tuple([testpt()...] for _ in 1:ntest)
# Generate the particles
test_points = florians_test_points();

function inside_torus(X,R₀=5.0,a=1.0)
    # Distance to axis
    dR = sqrt(X[1]^2 + X[2]^2) - R₀
    # Height above axis
    dZ = X[3]
    #If inside sign == -1, outside sign == 1, on surface == 0
    return sign(sqrt(dR^2 + dZ^2) - a)
end

exact_inside_outside = map(inside_torus, test_points);
computed_inside_outside = solid_angle_test(m,test_points);

@benchmark solid_angle_test($m,$test_points)




inout = solid_angle_test(m,test_points)


