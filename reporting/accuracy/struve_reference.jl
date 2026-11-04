#  (C) Copyright Nick Thompson 2026.
#  Use, modification and distribution are subject to the
#  Boost Software License, Version 1.0. (See accompanying file
#  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Regenerates the reference values in test/test_struve.cpp at the (nu, x) it lists:
#   julia struve_reference.jl ../../test/test_struve.cpp
# Needs SpecialFunctions (gamma of a BigFloat, which is MPFR's).
#
# H_nu(x) and L_nu(x) from the power series DLMF 11.2.1 and 11.2.2,
#   sum_k (-+1)^k (x/2)^(2k+nu+1) / (Gamma(k+3/2) Gamma(k+nu+3/2)),
# and their x-derivatives term by term, in BigFloat at the exact double arguments. The series converges for every x,
# but for H its terms reach about e^x / x before they cancel, so the working precision is chosen from the largest
# term and the computed value: 200 bits survive the cancellation. Each value is computed again with 128 more bits
# and the two must agree to 150 bits.

using SpecialFunctions
using Printf

function series(modified::Bool, nu::Float64, x::Float64)
    nub = BigFloat(nu)
    xb = BigFloat(x)
    x2 = xb / 2
    q = modified ? x2^2 : -x2^2
    ks = 0
    z = nu + 1.5
    if z <= 0 && isinteger(z)
        ks = Int(-z) + 1   # 1/Gamma(k + nu + 3/2) = 0 for k < ks
    end
    t = x2^(nub + 2ks + 1) / (gamma(BigFloat(ks) + BigFloat(3) / 2) * gamma(nub + ks + BigFloat(3) / 2))
    if !modified && isodd(ks)
        t = -t
    end
    s = zero(BigFloat)
    ds = zero(BigFloat)
    maxterm = zero(BigFloat)
    k = ks
    while true
        s += t
        ds += t * (2k + nub + 1) / xb
        maxterm = max(maxterm, abs(t))
        if k + nu + 1.5 > 0 && x^2 / 4 < (k + 1.5) * (k + nu + 1.5) / 4 && !iszero(s) &&
           abs(t) < abs(s) * BigFloat(2)^(-precision(BigFloat))
            break
        end
        t *= q / ((k + BigFloat(3) / 2) * (k + nub + BigFloat(3) / 2))
        k += 1
    end
    s, ds, maxterm
end

# f and f' with at least 150 correct bits.
function struve(modified::Bool, nu::Float64, x::Float64)
    p = 256
    while true
        s, ds, m = setprecision(BigFloat, p) do
            series(modified, nu, x)
        end
        lost = log2(m) - min(log2(abs(s)), log2(abs(ds * x)))
        if p < lost + 200
            p = ceil(Int, lost) + 264
            continue
        end
        s2, ds2, _ = setprecision(BigFloat, p + 128) do
            series(modified, nu, x)
        end
        @assert abs(s - s2) <= abs(s2) * BigFloat(2)^-150 "value not converged at nu = $nu, x = $x"
        @assert abs(ds - ds2) <= abs(ds2) * BigFloat(2)^-150 "derivative not converged at nu = $nu, x = $x"
        return s2, ds2
    end
end

fmt(v) = @sprintf("%.24e", v)
fmt50(v) = @sprintf("%.54e", v)

# The (nu, x) of the rows { nu, x, f, f' } of h_references and l_references, and of the struve_h/l(Real(nu), Real(x)) checks.
const src = read(ARGS[1], String)
table(name) = [(parse(Float64, m[1]), parse(Float64, m[2])) for m in
               eachmatch(r"\{ ([-0-9.]+), ([-0-9.]+), ", split(split(src, name)[2], "};")[1])]
const h_points = table("h_references[] = {")
const l_points = table("l_references[] = {")
const mp_points = [(m[1] == "l", parse(Float64, m[2]), parse(Float64, m[3])) for m in
                   eachmatch(r"struve_([hl])\(Real\(([-0-9.]+)\), Real\(([-0-9.]+)\), &d\), d,", src)]

function rows(modified, points)
    for (nu, x) in points
        f, df = setprecision(BigFloat, 256) do
            struve(modified, nu, x)
        end
        println("   { $(repr(nu)), $(repr(x)), $(fmt(f)), $(fmt(df)) },")
    end
end

println("// H_nu(x), H'_nu(x)")
rows(false, h_points)
println("// L_nu(x), L'_nu(x)")
rows(true, l_points)

println("// 55 digits, for cpp_bin_float_50")
for (modified, nu, x) in mp_points
    f, df = struve(modified, nu, x)
    println("   ", modified ? "L" : "H", "($(repr(nu)), $(repr(x))) = ", fmt50(f), "  ", fmt50(df))
end
