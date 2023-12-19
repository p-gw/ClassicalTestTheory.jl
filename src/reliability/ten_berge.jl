"""
    mu(m::AbstractMatrix, r::Int)
    mu(test::PsychometricTest, r::Int)
    mu(test::PsychometricTest, scale::Symbol, r::Int)

Calculate the lower bound of the reliability mu derived in $TENBERGE1978

## Notes
- If `r = 0` then mu is equivalent to Cronbach's alpha.
- If `r = 1` then mu is equivalent to Guttman's lambda₂.
"""
function mu(m::AbstractMatrix, r::Int)
    r >= 0 || throw(ArgumentError("r must be non-negative."))

    n = size(m, 2)
    C = cov(m)
    zerodiag!(C)

    st = var(sum(m, dims = 2))

    p_sum = zero(st)

    for h in Iterators.reverse(0:r)
        p_h = sum(c -> c^(2.0^h), C)

        if h == r
            p_h *= n / (n - 1)
        end

        if h == 0
            p_sum = p_sum + p_h
        else
            p_sum = sqrt(p_sum + p_h)
        end
    end

    return p_sum / st
end

@kwdef struct Mu <: ReliabilityMeasure
    r::Int
end

(method::Mu)(data) = mu(data, method.r)
name(r::Mu) = "Mu($(r.r))"

"""
    mu_up_to(r::Int)

Generate a vector of reliability measures `Mu` from 0 to `r`.

```jldoctest
julia> mu_up_to(2)
3-element Vector{Mu}:
 (::Mu) (generic function with 1 method)
 (::Mu) (generic function with 1 method)
 (::Mu) (generic function with 1 method)
```
"""
mu_up_to(r::Int) = [Mu(i) for i in 0:r]
