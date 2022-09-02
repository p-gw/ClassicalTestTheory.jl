"""
    find(t::AbstractTest, n; criterion=glb)

Perform an exhaustive search to find the optimal subtest of `t` with `n` items.

"""
function find(t, n; criterion=glb)
    is = vcat(trues(n), falses(nitems(t) - n))
    perms = multiset_permutations(is, length(is))
    n_perms_iter = Int(length(perms) / 2)

    optimal_subtest = is
    max_crit = -Inf
    @showprogress for perm in Iterators.take(perms, n_perms_iter)
        subtest, _ = split(t, perm)
        crit = criterion(subtest)
        if crit > max_crit
            max_crit = crit
            optimal_subtest = subtest
        end
    end
    return optimal_subtest
end
