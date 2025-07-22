using PointInPolygon
using StaticArrays

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




θₙ = 100
ζₙ = 500

verts,conns = triangulate_toroidal_surface(θₙ,ζₙ);

m = Mesh{Float64}(verts,conns);






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


inout_solidangle = solid_angle(m,test_points)

inout_winding = winding_number(m,test_points)



# @benchmark solid_angle($m,$test_points)
# @benchmark winding_number($m,$test_points)



# @profview solid_angle(m,test_points)