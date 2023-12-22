function zerodiag!(m)
    dims = size(m)
    foreach(i -> m[i, i] = zero(eltype(m)), 1:dims[2])
    return m
end

zerodiag(m) = m - diagm(diag(m))

function CTTPanel(args...; title = nothing)
    return Term.Panel(
        args...;
        title = isnothing(title) ? nothing : "{dim}" * title * "{/dim}",
        style = "dim",
        subtitle = "{dim}ClassicalTestTheory.jl{/dim}",
        subtitle_justify = :right,
        fit = true,
    )
end

const pretty_names = Dict(
    :lambda1 => "Guttman L₁",
    :lambda2 => "Guttman L₂",
    :lambda3 => "Guttman L₃ / Cronbachs α",
    :lambda4 => "Guttman L₄",
    :lambda5 => "Guttman L₅",
    :lambda6 => "Guttman L₆",
    :glb => "GLB",
)

prettify(s::Symbol) = pretty_names[s]


methodname(x) = typeof(x).name.name

quantilesfromlevel(l::Real) = ((1 - l) / 2, (1 + l) / 2)
prettyquantile(q::Real) = @sprintf("%.1f", q * 100) * "%-quantile"
