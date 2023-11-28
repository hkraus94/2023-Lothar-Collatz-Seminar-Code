function StefanTestCase(
    T_range,            #temperature range
    λ_range,            #jumping heat conductivity
    s_range,            #jumping volumetric heat conductivity
    T_phase_change,     #phase change teperature
    latent_heat         #latent heat
)
    T_cold, T_hot = T_range
    λ_cold, λ_hot = λ_range
    s_cold, s_hot = s_range

    α_hot  = λ_hot / s_hot
    α_cold = λ_cold / s_cold
    α_frac = α_hot / α_cold

    # the problem with and without latent heat is basically the same, differing by this factor
    lh_factor = !iszero(latent_heat) ? 1.0 / latent_heat : 1.0

    # Stefan numbers
    ST_hot  = s_hot * (T_hot - T_phase_change) * lh_factor
    ST_cold = s_cold * (T_phase_change - T_cold) * lh_factor

    nu = sqrt(α_frac)

    f(x) = ST_hot * nu * erfc(nu * x) * exp(nu^2 * x^2) -
        ST_cold * erf(x) * exp(x^2) -
        nu * erfc(nu * x) * erf(x) * sqrt(pi) * x *
        exp(nu^2 * x^2) * exp(x^2) * (!iszero(latent_heat))

    # nonlinear problem in different formulation to test if solution is good enough
    f2(x) = ST_hot / (exp(x^2) * erf(x)) - ST_cold / (nu * exp(nu^2 * x^2) * erfc(nu * x)) -
            x * sqrt(pi) * (!iszero(latent_heat))

    # the guess (0.7) should ideally depend on the model parameters
    x_left  = 0.0
    x_right = 0.0
    x_gain  = 0.1
    while f(x_left) * f(x_right) >= 0
        x_left   = x_right
        while isinf(f(x_right+x_gain)) || isnan(f(x_right+x_gain))
            x_gain *= 0.5
        end
        x_right += x_gain
    end
    λ = bisection(f, x_left, x_right)

    @show λ, f(λ), f2(λ)

    conds = [λ <= 0, abs(f(λ)) > 1e-10, abs(f2(λ)) > 1e-6]

    if any(conds)
        conds_violated = findall(conds)
        error("Conditions $conds_violated violated! NO FURTHER COMPUTATION!")
    end

    X(t) = 2λ * sqrt(α_hot * t)

    T_L(x, t) = T_hot - (T_hot - T_phase_change) * erf(x / 2sqrt(α_hot * t)) / erf(λ)
    T_R(x, t) = T_cold + (T_phase_change - T_cold) * erfc(x / 2sqrt(α_cold * t)) / erfc(λ * sqrt(α_frac))

    function T(x::Number, t)
        if iszero(t)
            return x <= X(t) ? T_hot : T_cold
        else
            return x <= X(t) ? T_L(x, t) : T_R(x, t)
        end
    end

    return (
        solution = T,
        interface = X
    )
end