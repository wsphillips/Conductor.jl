
abstract type Geometry end

struct Sphere <: Geometry
    radius::Any
end

function Sphere(; radius)
    Sphere(radius)
end

struct Point <: Geometry end

struct Unitless <: Geometry
    value::Any
end

struct Cylinder <: Geometry
    radius::Any
    height::Any
    open_ends::Bool
end

function Cylinder(; radius, height, open_ends = true)
    return Cylinder(radius, height, open_ends)
end

height(x::Geometry) = isdefined(x, :height) ? getfield(x, :height) : nothing
radius(x::Geometry) = getfield(x, :radius)
radius(::Union{Point, Unitless}) = 0.0

area(x::Sphere) = ustrip(Float64, cm^2, 4 * π * radius(x)^2)
function area(x::Cylinder)
    ustrip(Float64, cm^2, 2 * π * radius(x) * (height(x) + (x.open_ends ? 0µm : radius(x))))
end
area(::Point) = 1.0
area(x::Unitless) = x.value
