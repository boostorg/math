#  (C) Copyright Nick Thompson 2026.
#  Use, modification and distribution are subject to the
#  Boost Software License, Version 1.0. (See accompanying file
#  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Generates the reference table in test/test_spheroidal.cpp:
#   julia spheroidal_reference.jl
# Needs Nemo (Arb) for the spherical Bessel functions; everything else is MPFR (BigFloat).
#
# Eigenvalues: Sturm bisection on the symmetrized tridiagonal matrix of A&S 21.7.3.
# Coefficients d_k: ratio recurrences from both ends, d_{n-m} = 1, then the Meixner-Schäfke
# normalization and the sign of DLMF 30.4.
# Angular functions: Σ d_k P_{m+k}^m(η) (Ferrers functions by recurrence).
# Radial functions: the expansions in j_{m+k}(cξ) and y_{m+k}(cξ), normalized by Σ d_k (2m+k)!/k!.
# At 384 bits their cancellation is harmless. Below ξ = 2 the Neumann series converges slowly, so
# R2 (and the oblate R1 at ξ = 0) is continued from ξ = 2 by a Taylor-series method; every row
# is checked against the Wronskian to 40 digits.

using Nemo

# A&S 21.7.3, with s = +c^2 (prolate) or -c^2 (oblate)
alpha(m, k, s) = (2m + k + 2) * (2m + k + 1) * s / ((2m + 2k + 3) * (2m + 2k + 5))
beta(m, k, s)  = (m + k) * (m + k + 1) + (2 * (m + k) * (m + k + 1) - 2m^2 - 1) * s / ((2m + 2k - 1) * (2m + 2k + 3))
gamma_(m, k, s) = k * (k - 1) * s / ((2m + 2k - 3) * (2m + 2k - 1))

# Number of eigenvalues of the symmetrised tridiagonal block below x (Sturm count).
function sturm_count(diag, offsq, x)
    cnt = 0
    q = diag[1] - x
    q < 0 && (cnt += 1)
    for i in 2:length(diag)
        q == 0 && (q = eps(x) * (abs(x) + 1))
        q = diag[i] - x - offsq[i-1] / q
        q < 0 && (cnt += 1)
    end
    cnt
end

function eigenvalue(m, n, s, N)
    p = (n - m) % 2
    ks = [p + 2j for j in 0:N-1]
    diag = [beta(m, k, s) for k in ks]
    offsq = [alpha(m, ks[i], s) * gamma_(m, ks[i+1], s) for i in 1:N-1]
    j = (n - m - p) ÷ 2
    r = maximum(abs, diag) + 2 * maximum(sqrt ∘ abs, offsq)
    lo, hi = -r, r
    for _ in 1:(precision(BigFloat) + 64)
        mid = (lo + hi) / 2
        if sturm_count(diag, offsq, mid) > j
            hi = mid
        else
            lo = mid
        end
    end
    (lo + hi) / 2
end

# Flammer d_k for k ≡ n-m (mod 2), k = p, p+2, ..., p+2(N-1); d_{n-m} = 1 before normalisation.
function dcoeffs(m, n, s, λ, N)
    p = (n - m) % 2
    ks = [p + 2j for j in 0:N-1]
    d = Dict{Int,BigFloat}()
    i0 = (n - m - p) ÷ 2 + 1
    d[ks[i0]] = one(BigFloat)
    # backward ratios ρ_k = d_k / d_{k-2}, minimal solution towards large k
    ρ = zero(BigFloat)
    ratios = Dict{Int,BigFloat}()
    for i in N:-1:i0+1
        k = ks[i]
        ρ = -gamma_(m, k, s) / (beta(m, k, s) - λ + alpha(m, k, s) * ρ)
        ratios[k] = ρ
    end
    for i in i0+1:N
        d[ks[i]] = d[ks[i-1]] * ratios[ks[i]]
    end
    # forward ratios σ_k = d_k / d_{k+2}
    σ = zero(BigFloat)
    sig = Dict{Int,BigFloat}()
    for i in 1:i0-1
        k = ks[i]
        σ = -alpha(m, k, s) / (beta(m, k, s) - λ + gamma_(m, k, s) * σ)
        sig[k] = σ
    end
    for i in i0-1:-1:1
        d[ks[i]] = d[ks[i+1]] * sig[ks[i]]
    end
    ks, [d[k] for k in ks]
end

# Associated Legendre P_l^m(0) (even l-m) or dP_l^m/dη(0) (odd l-m), no Condon-Shortley phase.
function legendre_at_zero(l, m)
    if (l - m) % 2 == 0
        (-1)^((l - m) ÷ 2) * factorial(big(l + m)) / (big(2)^l * factorial(big((l - m) ÷ 2)) * factorial(big((l + m) ÷ 2)))
    else
        (-1)^((l - m - 1) ÷ 2) * factorial(big(l + m + 1)) / (big(2)^l * factorial(big((l - m - 1) ÷ 2)) * factorial(big((l + m + 1) ÷ 2)))
    end
end


# P_l^m(η) and derivative, l = m .. m+L, no Condon-Shortley phase.
function legendre_table(m, L, η)
    P = Vector{BigFloat}(undef, L + 1)
    w = 1 - η^2
    pmm = one(BigFloat)
    for i in 1:m
        pmm *= (2i - 1)
    end
    pmm *= sqrt(w)^m
    P[1] = pmm
    L >= 1 && (P[2] = η * (2m + 1) * pmm)
    for i in 2:L
        l = m + i - 1
        P[i+1] = ((2l + 1) * η * P[i] - (l + m) * P[i-1]) / (l - m + 1)
    end
    dP = similar(P)
    for i in 1:L+1
        l = m + i - 1
        prev = i == 1 ? zero(BigFloat) : P[i-1]
        dP[i] = ((l + m) * prev - l * η * P[i]) / w
    end
    P, dP
end


# Spherical Bessel j_l / y_l for l = 0..L at real z, via Arb.
function sph_bessel(L, z::BigFloat, kind)
    prec = precision(BigFloat) + 64
    C = AcbField(prec)
    R = ArbField(prec)
    zc = C(R(z))
    f = kind == 1 ? bessel_j : bessel_y
    pre = sqrt(C(R(pi)) / (2 * zc))
    [BigFloat(midpoint(real(pre * f(C(l) + C(1) // C(2), zc)))) for l in 0:L]
end

function radial(m, n, s, ξ, N)
    λ = eigenvalue(m, n, s, N)
    ks, d = dcoeffs(m, n, s, λ, N)
    σ = s > 0 ? 1 : -1
    c = sqrt(abs(s))
    z = c * ξ
    L = m + maximum(ks) + 1
    j = sph_bessel(L, z, 1)
    y = sph_bessel(L, z, 2)
    w = [d[i] * factorial(big(2m + ks[i])) / factorial(big(ks[i])) for i in eachindex(ks)]
    sgn = [iseven((ks[i] - n + m) ÷ 2) ? 1 : -1 for i in eachindex(ks)]
    norm = sum(w)
    f = ((ξ^2 - σ) / ξ^2)^(BigFloat(m) / 2)
    dlogf = σ * m / (ξ * (ξ^2 - σ))
    out = BigFloat[]
    for b in (j, y)
        S = sum(sgn[i] * w[i] * b[m+ks[i]+1] for i in eachindex(ks))
        # d/dz b_l = b_{l-1} - (l+1)/z b_l  (b_{-1} not needed: l = m+k >= 0, use b_1 relation for l=0)
        dS = sum(begin
                     l = m + ks[i]
                     db = l == 0 ? -b[2] : b[l] - (l + 1) / z * b[l+1]
                     sgn[i] * w[i] * db
                 end for i in eachindex(ks))
        R = f * S / norm
        Rd = (dlogf * f * S + f * c * dS) / norm
        push!(out, R, Rd)
    end
    λ, out...
end

# ---- Taylor-series continuation of the radial ODE ----
# (ξ²-1)² R'' + 2ξ(ξ²-1) R' - [(λ - sξ²)(ξ²-1) + m²] R = 0   (s = c² prolate)
polymul(a, b) = [sum(a[i] * b[k-i+2] for i in max(1, k - length(b) + 2):min(length(a), k + 1)) for k in 0:length(a)+length(b)-2]
function ode_coeffs(a, m, s, λ)
    σ = s > 0 ? 1 : -1
    ξ = [a, one(a)]                       # a + t
    q = polymul(ξ, ξ); q[1] -= σ          # (a+t)² - σ
    P2 = polymul(q, q)
    P1 = 2 .* polymul(ξ, q)
    lin = -abs(s) .* polymul(ξ, ξ); lin[1] += λ
    P0 = -polymul(lin, q); P0[1] -= σ * m^2
    P2, P1, P0
end
# Step from a to b (|b-a| < a-1), returning (R(b), R'(b)).
function taylor_step(a, b, R, Rd, m, s, λ)
    P2, P1, P0 = ode_coeffs(a, m, s, λ)
    h = b - a
    r = [R, Rd]
    tol = eps(BigFloat)^2
    k = 0
    while true
        # coefficient of t^k
        acc = zero(BigFloat)
        for j in 0:4
            if j > 0 && k - j + 2 >= 0 && j + 1 <= length(P2)
                acc += P2[j+1] * (k - j + 2) * (k - j + 1) * r[k-j+3]
            end
            if k - j + 1 >= 0 && j + 1 <= length(P1)
                acc += P1[j+1] * (k - j + 1) * r[k-j+2]
            end
            if k - j >= 0 && j + 1 <= length(P0)
                acc += P0[j+1] * r[k-j+1]
            end
        end
        push!(r, -acc / (P2[1] * (k + 2) * (k + 1)))
        k += 1
        if k > 20 && abs(r[end] * h^(k + 1)) + abs(r[end-1] * h^k) < tol * (abs(R) + abs(Rd) + 1)
            break
        end
        k > 20000 && error("Taylor series did not converge")
    end
    Rb = sum(r[i] * h^(i - 1) for i in eachindex(r))
    Rdb = sum((i - 1) * r[i] * h^(i - 2) for i in 2:length(r))
    Rb, Rdb
end
function continue_radial(ξ0, ξ, R, Rd, m, s, λ)
    a = ξ0
    while a != ξ
        h = s > 0 ? (a - 1) / 2 : sqrt(a^2 + 1) / 2
        h = min(h, 1 / (2 * sqrt(abs(s))) + h / 4)
        b = ξ < a ? max(ξ, a - h) : min(ξ, a + h)
        R, Rd = taylor_step(a, b, R, Rd, m, s, λ)
        a = b
    end
    R, Rd
end



# DLMF 30.4: ||Ps|| = ||P_n^m||, sign of Ps(0) [even] or Ps'(0) [odd] equals that of P_n^m.
function ms_normalise!(d, ks, m, n)
    t = sum(d[i] * legendre_at_zero(m + ks[i], m) for i in eachindex(ks))
    nrm = sum(d[i]^2 * 2 / big(2m + 2ks[i] + 1) * factorial(big(2m + ks[i])) / factorial(big(ks[i])) for i in eachindex(ks))
    target = 2 / big(2n + 1) * factorial(big(n + m)) / factorial(big(n - m))
    d .*= sqrt(target / nrm) * sign(t / legendre_at_zero(n, m))
    d
end
function angular_ms(m, n, s, η, N)
    λ = eigenvalue(m, n, s, N)
    ks, d = dcoeffs(m, n, s, λ, N)
    ms_normalise!(d, ks, m, n)
    P, dP = legendre_table(m, maximum(ks) + 1, η)
    cs = isodd(m) ? -1 : 1           # Condon-Shortley phase
    λ, cs * sum(d[i] * P[ks[i]+1] for i in eachindex(ks)), cs * sum(d[i] * dP[ks[i]+1] for i in eachindex(ks))
end

# Reference rows for test/test_spheroidal.cpp: Arb (Nemo) + MPFR, see swf_ref.jl.
setprecision(BigFloat, 384)
using Printf
fmt(x) = x == 0 ? "0" : @sprintf("%.24e", x)
eig = Set{Tuple{String,Int,Int,String}}()
rows = String[]
function row(fam, kind, n, m, cstr, xstr)
    c = parse(BigFloat, cstr); x = parse(BigFloat, xstr)
    s = fam == "prolate" ? c^2 : -c^2
    N = (n - m) ÷ 2 + Int(ceil(Float64(c))) + 150
    push!(eig, (fam, n, m, cstr))
    if kind == "angular"
        λ, S, dS = angular_ms(m, n, s, x, N)
        d2 = (2x * dS - (λ - s * x^2 - m^2 / (1 - x^2)) * S) / (1 - x^2)
        vals = (S, dS, d2)
    else
        λ = eigenvalue(m, n, s, N)
        σ = s > 0 ? 1 : -1
        _, A1, A1d, A2, A2d = radial(m, n, s, big(2), N)
        if x >= 2
            _, R1, R1d, R2, R2d = radial(m, n, s, x, N)
        else
            R2, R2d = continue_radial(big(2), x, A2, A2d, m, s, λ)
            if x == 0 || (σ > 0 && x == 1)
                R1, R1d = σ < 0 ? continue_radial(big(2), x, A1, A1d, m, s, λ) : (radial(m, n, s, x + big(10)^-60, N)[2:3])
            else
                _, R1, R1d, _, _ = radial(m, n, s, x, N)
            end
        end
        q = x^2 - σ
        W = (R1 * R2d - R1d * R2) * c * q
        abs(W - 1) > big(10)^-40 && error("Wronskian check failed $fam $n $m $cstr $xstr")
        dd(R, Rd) = -(2x * q * Rd + ((c^2 * x^2 - λ) * q - σ * m^2) * R) / q^2
        vals = kind == "radial_1" ? (R1, R1d, dd(R1, R1d)) : (R2, R2d, dd(R2, R2d))
    end
    push!(rows, "   { $(fam == "oblate"), quantity::$kind, $n, $m, $cstr, $xstr, $(fmt(vals[1])), $(fmt(vals[2])), $(fmt(vals[3])) },")
end
for (n, m, c, x) in [(2,1,"1","0.3"), (7,3,"2","0.5"), (5,2,"8","-0.7"), (3,0,"20","0.95"), (4,1,"50","0.999"), (0,0,"50","0.3")]
    row("prolate", "angular", n, m, c, x)
end
for (n, m, c, x) in [(0,0,"1","1.5"), (2,1,"2","1.5"), (1,0,"2","1.00000001"), (3,1,"8","1.001"), (5,2,"20","1.3"), (0,0,"50","1.05"), (15,5,"20","20"), (4,0,"50","5")]
    row("prolate", "radial_1", n, m, c, x); row("prolate", "radial_2", n, m, c, x)
end
for (n, m, c, x) in [(0,0,"50","0"), (5,5,"50","0.3"), (6,5,"20","-0.7"), (2,1,"2","0.95"), (3,0,"8","0.2")]
    row("oblate", "angular", n, m, c, x)
end
for (n, m, c, x) in [(0,0,"1","1.5"), (2,1,"2","0.3"), (2,0,"0.7","0.95"), (5,5,"50","0.01"), (7,5,"50","0.2"), (15,5,"0.5","3"), (3,2,"8","0"), (0,0,"20","1")]
    row("oblate", "radial_1", n, m, c, x); row("oblate", "radial_2", n, m, c, x)
end
println("// { oblate, quantity, n, m, c, x, f(x), f'(x), f''(x) }")
foreach(println, rows)
println("// eigenvalues: { oblate, n, m, c, DLMF lambda }")
for (fam, n, m, cstr) in sort(collect(eig))
    c = parse(BigFloat, cstr); s = fam == "prolate" ? c^2 : -c^2
    λ = eigenvalue(m, n, s, (n - m) ÷ 2 + Int(ceil(Float64(c))) + 150)
    println("   { $(fam == "oblate"), $n, $m, $cstr, $(fmt(λ - s)) },")
end
