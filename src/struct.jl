struct BoxCoxTransformation{T <: AbstractFloat, N <: Real}
    value::T
    details::UnivariateOptimizationResults
    method::Symbol
    dat::Vector{N}
end

Base.propertynames(::BoxCoxTransformation) =
            (:lambda, :dat, :method, :details, :loglik, :transformed)

function Base.getproperty(bc::BoxCoxTransformation, s::Symbol)
    if s === :lambda || s === :λ
        return bc.value
    elseif s === :𝐱
        return bc.dat
    elseif s === :transformed
        return transform(bc.𝐱, bc.λ)
    elseif s === :loglik
        return log_likelihood(bc.dat, bc.λ; method = bc.method)
    else
        return getfield(bc, s)
    end
end

log_likelihood(bc::BoxCoxTransformation) = log_likelihood(bc.dat, bc.λ; method = bc.method)
log_likelihood(bc::BoxCoxTransformation, λ) = log_likelihood(bc.dat, λ; method = bc.method)

transform(bc::BoxCoxTransformation; kwargs...) = transform(bc.𝐱, bc.λ; kwargs...)
transform(bc::BoxCoxTransformation, λ::Real; kwargs...) = transform(bc.𝐱, λ; kwargs...)

"""
    BoxCoxTransformation(𝐱::Vector{<:Real}; interval = (-2.0, 2.0),
        method = :geomean)

Create a BoxCoxTransformation. The method wraps `BoxCoxTrans.lambda`` but
returns a BoxCoxTransformation

For details see: [`lambda`](@ref)
"""
function BoxCoxTransformation(𝐱::Vector{<:Real};
    interval = (-2.0, 2.0), method = :geomean, kwargs...)
    bc = lambda(𝐱; interval, method, kwargs...)
    return BoxCoxTransformation(bc.value, bc.details, method, 𝐱)
end

"""
    confint(bc::BoxCoxTransformation; alpha=0.05)

Return confidence interval for the estimated maximum power parameter λ of a
given BoxCoxTransformation.

Note: The confidence interval is the range of λs for which the
log-likelihood(λ) >= log_likelihood(λ_max) - 0.5 * quantile(Chisq(1), 1-alpha)

Reference: Linear Models With R by Julian Faraway Section 8.1.
"""
function confint(bc::BoxCoxTransformation; alpha = 0.05)
    interval = (bc.details.initial_lower, bc.details.initial_upper)
    τ = log_likelihood(bc, bc.λ) - 0.5 * quantile(Chisq(1), 1 - alpha) # target log-likelihood
    dist_fnc(x) = abs(log_likelihood(bc, x) - τ)
    res_upper = optimize(dist_fnc, bc.λ, interval[2])
    res_lower = optimize(dist_fnc, interval[1], bc.λ)
    return (minimizer(res_lower), minimizer(res_upper))
end


function Base.show(io::IO, ::MIME"text/plain", x::BoxCoxTransformation)
	println(io, "BoxCoxTransformation: λ = " * string(x.value))
	#println(io, x.details)
end;