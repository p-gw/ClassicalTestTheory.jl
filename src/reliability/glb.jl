"""
    glb(m::AbstractMatrix)

Return the greatest lower bound estimate (glb) of the reliability as described in
$WOODHOUSE1977
"""
function glb(m::AbstractMatrix)
    n = size(m, 2)

    C = cov(m)
    C̃ = zerodiag(C)
    upr = diag(C)
    lwr = zeros(n)

    model = JuMP.Model(SCS.Optimizer)

    set_silent(model)
    set_string_names_on_creation(model, false)

    @variable(model, y[1:n])
    @expression(model, A, Symmetric(C̃ + diagm(y)))

    @objective(model, Min, sum(y))
    @constraint(model, lwr .<= y .<= upr)
    @constraint(model, A in PSDCone())

    optimize!(model)

    if termination_status(model) == OPTIMAL
        sum_y = sum(value.(y))
        return (sum(C̃) + sum_y) / sum(C)
    else
        error("something went wrong")  # TODO: Fix error message
    end
end

struct GLB <: ReliabilityMeasure end
(method::GLB)(data) = glb(data)
