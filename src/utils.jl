function zerodiag!(m)
    dims = size(m)
    foreach(i -> m[i, i] = zero(eltype(m)), 1:dims[2])
    return m
end

zerodiag(m) = m - diagm(diag(m))
