################################################################################
## The code is based on the power spectrum analysis due to variation in width ##
## and spacings done by A. S. Dogra. Written by @Pritipriya_dasbehera         ##
################################################################################


using Plots, FFTW
using Base: mod2pi

function widths(om::Float64, w::Float64, delta; tol=0.2)
    base = exp(-w^2 * om^2 / (1 + 2 * (delta^2) * (om^2))) / sqrt(1 + 2 * (delta^2) * (om^2))

    # distance from nearest multiple of 2π
    dist = abs(mod2pi(om))
    dist = min(dist, 2π - dist)  # handle wrap-around

    # Gaussian mask centered at multiples of 2π
    mask = exp(-(dist^2) / (2 * tol^2))

    return base * mask
end


# --- Single omega ---
function widths(om::Vector{Float64}, w::Float64, delta; tol=0.2)
    res = zeros(length(om))
    for i in eachindex(om)
        res[i] = widths(om[i], w, delta; tol=tol)
    end
    return res
end

function test_widths_plot()
    w = 0.1
    delta = 0
    tol = 0.01

    # sample range over several multiples of 2π
    om = collect(-10π:0.01:10π)
    y = widths(om, w, delta; tol=tol)

    # plot
    plot(om, y,
         lw=2,
         label="δ = $delta, ω = $om, tol = $tol",
         xlabel="w",
         ylabel="widths(w)",
         title="Widths with δ-comb Mask",
         legend=:topright,
         framestyle=:box)

    # mark multiples of 2π for visual reference
    for n in -5:5
        vline!([2π * n], color=:red, lw=0.5, ls=:dot, label=false)
    end

    display(current())
end


function plt()
#     om = 2pi.*collect(-10:1:10)
    w = 0.16
    deltas = [0.0, 0.02, 0.05]
    x = collect(-500pi:pi/8:500pi)

    xmin, xmax = -25, 25
    idf = findall(f -> xmin ≤ f ≤ xmax, x)

    dx = x[2] - x[1]
    n = length(x)

    # Frequency axis for FFT (centered)
    freqs = fftshift(fftfreq(n, 1/dx)) *2pi

    # Choose frequency window
    fmin, fmax = 0, 2
    idx = findall(f -> fmin ≤ f ≤ fmax, freqs)

    # Set up nice color palette and layout
    palt = palette(:viridis, length(deltas)+3)
    alpha = [0.7, 0.5, 0.4]
    plt = plot(layout = (2, 1),
               size=(600, 500),
#                 xlabel="f", ylabel="magnitude of f",
               title="Widths and their Fourier Transforms",
               bglegend=RGBA(1,1,1,0.6),
               legend=:topright, lw=1, framestyle=:box)

    # Time-domain plots
    for (i, δ) in enumerate(deltas)
        y = widths(x, w, δ;tol= 0.01)
        plot!(plt[1], x[idf], y[idf], color=palt[i], alpha = alpha[i],
#               marker=:+,
            label="δ = $(δ)", title="(a) Frequency domain", )

        # Compute normalized FFT and shift
        Y = fftshift(abs.(fft(y)))
        Y ./= maximum(Y)
        plot!(plt[2], freqs[idx], Y[idx], alpha = alpha[i],
#               marker=:+,
color=palt[i], label="δ = $(δ)", title=" (b) Real space")
    end

    # Beautify both subplots
    plot!(plt[1], xlabel="Frequency (ω)", ylabel="Magnitude", legend=:topright)
    plot!(plt[2], xlabel="Position (x)", ylabel="Correlation", legend=:topright)
#     plot!(plt, )
#     println(length(idx))
#     println(length(freqs))
    display(plt)
#     savefig(plt, "var-width.svg")
savefig(plt, "var-width-temp.png")
end

function normal(om, mu, sigma)
    σ = (sigma == 0) ? 1e-10 : sigma
    return exp( -((om - mu)^2) / (2 * σ^2) )
end

function spacings(om::Float64, w::Float64, delta; K=5)
    res = 0
    for k in -K:1:K
        res += normal(om, 2pi*k, 2pi*abs(k)*delta) * exp( - w^2 * om^2)
    end
    return res
end

function spacings(om::Vector{Float64}, w::Float64, delta; K=20)
    res = zeros(length(om))
    for i in eachindex(om)
        res[i] += spacings(om[i], w, delta; K=K)
    end
    return res
end

# --- Plotting function ---
function plt_spectrum()
    w = 0.16               # Gaussian window
    deltas = [0, 0.02, 0.05]   # sampling uncertainty
    x = collect(-100pi:pi/32:100pi)  # "time" axis

    xmin, xmax = -25, 25
    idf = findall(f -> xmin ≤ f ≤ xmax, x)

    dx = x[2] - x[1]
    n = length(x)

    # Frequency axis for FFT (centered)
    freqs = fftshift(fftfreq(n, 1/dx)) * 2π

    # Frequency window for plotting
    fmin, fmax = 0, 4
    idx = findall(f -> fmin ≤ f ≤ fmax, freqs)

    # Color palette and alpha
    palt = palette(:viridis, length(deltas)+1)
    alpha = [0.7, 0.6, 0.5]

    plt = plot(layout=(2,1), size=(600,500),
               title="Spectrum vs Fourier Transform",
               bglegend=RGBA(1,1,1,0.6),
               legend=:topright, lw=2, framestyle=:box)

    for (i, delta) in enumerate(deltas)
        # Time-domain "spectrum" (function of ω)
        y = spacings(x, w, delta)

        plot!(plt[1], x[idf], y[idf], color=palt[i], alpha=alpha[i],
              label="δ = $(delta)", title="(a) Frequency domain", )

        # FFT of y (centered and normalized)
        Y = fftshift(abs.(fft(y)))
        Y ./= maximum(Y)
        plot!(plt[2], freqs[idx], Y[idx], color=palt[i], alpha=alpha[i],
              label="δ = $(delta)", title=" (b) Real space")
    end


    plot!(plt[1], xlabel="Frequency (ω)", ylabel="Magnitude", legend=:topleft)
    plot!(plt[2], xlabel="Position (x)", ylabel="Correlation", legend=:topleft)

    display(plt)
#     savefig(plt, "var-spacing.svg")
    savefig(plt, "var-spacing-temp.png")
end

# --- Run the test ---
# test_gaussian_fft()
