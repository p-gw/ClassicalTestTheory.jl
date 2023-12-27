using PrecompileTools

@setup_workload begin
    m = rand(0:100, 20, 4)

    @compile_workload begin
        # reliability measures
        methods = [GLB(), GUTTMAN_METHODS..., KR20(), KR21(), mu_up_to(1)...]
        reliability(m, methods)

        result = reliability(m, Alpha())

        data(result)
        sampling_method(result)
        ci_method(result)
        estimate(result)
        estimate(result, "Alpha")
        confint(result)
        confint(result, "Alpha")
        bias(result)
        bias(result, "Alpha")
        stderror(result)
        stderror(result, "Alpha")
        samples(result, "Alpha")

        # item statistics
        ia = itemanalysis(m)

        # find
        for method in methods
            find(m, 2, method)
        end
    end
end
