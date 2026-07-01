
include("ShearHeatingShearBands.jl")

Diff    =   [
    :explicit,
    :implicit,
    :dc,
]

θ   =   [
    0,
    0.5,
    1.0
]

Adv     =   [
    :upwind,
    :slf,
    :semilag,
]

style   =   [
    :moving,
    :fixed,
]

for i in eachindex(θ)
    @show θ[i] 
    for j in eachindex(Adv)
        @show Adv[j]
        ShearHeatingShearBands(Diff[3],θ[i],Adv[j],style[1])
    end
end

