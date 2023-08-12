#These four convenience functions extract the state variable from the state vector
#It is assumed the layout of the vector u is u = [v_x, v_y, x, y]
state_v_x(u::Vector{Float64}) = u[1]
state_v_y(u::Vector{Float64}) = u[2]
state_x(u::Vector{Float64}) = u[3]
state_y(u::Vector{Float64}) = u[4]

"""
Computes the RHS for the one body problem. 
"""
function df_dt_one_body(u::Vector{Float64}, t::Float64)::Vector{Float64}
    M, G = 1, 1 #We take these constants as normalized. Naturally they would need to be set for physical values.
    r = sqrt(state_x(u)^2 + state_y(u)^2)
    return [-M*G*state_x(u)/r^3, -M*G*state_y(u)/r^3, state_v_x(u), state_v_y(u)]
end;

using Plots, Measures

function plot_solution( t::AbstractArray{T}, 
                        u::Vector{Vector{Float64}}; 
                        title::String = "",
                        label::Union{String, Bool} = false) where T
    x, y, v_x, v_y = state_x.(u), state_y.(u), state_v_x.(u), state_v_y.(u)
 
    #"Energy"
    r = @. sqrt(x^2 + y^2)
    E = @. 0.5*(v_x^2 + v_y^2) - 1.0/r

    p1 = plot(  x, y, label = label, xlabel= "X", ylabel = "Y",
                title = title*" (position)", aspectratio=1,legend=:topleft,ylim=(-7,7))
    scatter!([0], [0], ms=15, msw=0, c=:orange, shape =:star, label="Sun")
    scatter!([x[1]], [y[1]], ms=5, msw=0, c=:blue, shape =:circle, label="Earth initial position")

    p4 = plot(  t, E, xlabel = "Time", ylabel = "Energy",
                label = label, title = title*" (energy)")
    plot(p1, p4, margin = 10mm,size=(800,400))
end;

# Euler's method implementation
function euler_step(u::Vector{Float64}, h::Float64)::Vector{Float64}
    du = df_dt_one_body(u, 0.0)  # Calculate the derivatives at t = 0.0 (not used here)
    return u + h * du
end

# Initial conditions
u0 = [0.0, 2.0, 1.0, 0.0]  # [v_x, v_y, x, y]

# Time range
t_start = 0.0
t_end = 10.0
h_values = [0.01, 0.001, 0.0001]

# Arrays to store solutions
solutions = Dict()

for h in h_values
    num_steps = Int((t_end - t_start) / h) + 1
    t_values = LinRange(t_start, t_end, num_steps)
    u_values = [u0]
    
    u_current = u0
    for i in 1:num_steps
        u_next = euler_step(u_current, h)
        push!(u_values, u_next)
        u_current = u_next
    end
    
    solutions[h] = (t_values, u_values)
end

# Plotting the solutions
for h in h_values
    t_values, u_values = solutions[h]
    title = "Euler Method (h = $h)"
    
    display(plot_solution(t_values, u_values, title = title, label = "h = $h"))
end
