using Plots
using FFTW
using LinearAlgebra
using Distributions
using Interpolations
using ProgressMeter
using Printf

######################################################################################################################
# Define parameters

N_X = 101 # number of points in x direction
N_Y = 301 # number of points in y direction
Lx = 1.0 # length of grid in x direction
Ly = 3.0 # length of grid in y direction
μ = 0.0001 # viscosity
Δt = 0.0001 # timestep
N_TIME_STEPS = 900 # number of simulation timesteps
FORCE_STRENGTH = 200 # multiplicative force factor

######################################################################################################################
# Define constant variables
Δx = Lx / (N_X-1)
Δy = Ly / (N_Y-1)

x_values = 0.0:Δx:Lx
y_values = 0.0:Δy:Ly

x_coordinates = repeat(x_values, 1, N_Y)
y_coordinates = repeat(y_values', N_X, 1)

k_vals_x = fftfreq(N_X) .* N_X
k_vals_y = fftfreq(N_Y) .* N_Y

k_x = repeat(k_vals_x, 1, length(k_vals_y))
k_y = repeat(k_vals_y', length(k_vals_x), 1)

k_norm = [norm([k_x, k_y]) for k_x in k_vals_x, k_y in k_vals_y]
k_norm[iszero.(k_norm)] .= 1.0

k_x ./= k_norm
k_y ./= k_norm

filter = exp.(- Δt .* μ .* k_norm.^2)

force_x0 = FORCE_STRENGTH * (
    # - exp.( @. (-(x_coordinates .- 0.5).^2 - (y_coordinates .- 0.5).^2) / 0.005) 
    + exp.( @. (-(x_coordinates .- 0.5).^2 - (y_coordinates .- 0.5).^2) / 0.005)
)
force_y0 = FORCE_STRENGTH * (
    + exp.( @. (-(x_coordinates .- 0.5).^2 - (y_coordinates .- 0.5).^2) / 0.001)
    # - exp.( @. (-(x_coordinates .- 0.3).^2 - (y_coordinates .- 0.7).^2) / 0.001)
)

######################################################################################################################

function interpolate(Z₀::Matrix{Float64}, X::Matrix{Float64}, Y::Matrix{Float64})::Matrix{Float64}
    itp = LinearInterpolation((x_values, y_values), Z₀)
    interpolated = itp.(X, Y)
    return interpolated
end 

function past_position(R₀::Matrix{Float64}, V::Matrix{Float64}, ΔL)::Matrix{Float64}
    return mod1.(R₀ - Δt * V, ΔL)
end

function main()

    # preallocate variables
    force_x::Matrix{Float64} = zeros(Float64, N_X, N_Y)
    force_y::Matrix{Float64} = zeros(Float64, N_X, N_Y)

    velocity_x::Matrix{Float64} = zeros(Float64, N_X, N_Y)
    velocity_y::Matrix{Float64} = zeros(Float64, N_X, N_Y)

    advection_x::Matrix{Float64} = zeros(Float64, N_X, N_Y)
    advection_y::Matrix{Float64} = zeros(Float64, N_X, N_Y)

    freqs_x::Matrix{Complex{Float64}} = zeros(Complex{Float64}, N_X, N_Y)
    freqs_y::Matrix{Complex{Float64}} = zeros(Complex{Float64}, N_X, N_Y)
    projection::Matrix{Complex{Float64}} = zeros(Complex{Float64}, N_X, N_Y)

    @showprogress for iter = 1:N_TIME_STEPS
        
        # give the force a direction
        force_x = sin(0*pi + 0.05*sin(iter*0.1)) * force_x0
        force_y = cos(0*pi + 0.05*cos(iter*0.1)) * force_y0

        # apply the force
        velocity_x +=  force_x
        velocity_y +=  force_y

        # advection
        advection_x = past_position(x_coordinates, velocity_x, Lx)
        advection_y = past_position(y_coordinates, velocity_y, Ly)

        velocity_x = interpolate(velocity_x, advection_x, advection_y)
        velocity_y = interpolate(velocity_y, advection_x, advection_y)

        # filter high frequencies
        freqs_x = fft(velocity_x)
        freqs_y = fft(velocity_y)

        freqs_x .*= filter
        freqs_y .*= filter

        # make the fluid mass conserving (incompressible): 
        # project velocity to k direction (Helmholtz decomposition)
        # k vectors are already normalized
        projection = freqs_x .* k_x + freqs_y .* k_y
        freqs_x -= projection .* k_x
        freqs_y -= projection .* k_y

        # undo fourier transform
        velocity_x = real(ifft(freqs_x))
        velocity_y = real(ifft(freqs_y))

        # subtract mean velocity
        velocity_x .-= mean(vec(velocity_x))
        velocity_y .-= mean(vec(velocity_y))

        # compute vorticity to plot
        du_dy = diff(velocity_x, dims=2)[2:end, :]
        dv_dx = diff(velocity_y, dims=1)[:, 2:end]
        vorticity = du_dy - dv_dx

        # plot frame results
        theme(:dark)
        fig = heatmap(
            # vorticity,
            velocity_x.^2 + velocity_y.^2,
            c=:seaborn_icefire_gradient, 
            # clim=(-100,100), 
            # aspect_ratio=Ly/Lx, 
            size=(Ly*500,Lx*500),
            ticks=false
            )
        display(fig)
        # savefig(fig, "figures//$(@sprintf("%03d", iter)).png")

    end

end

main()
