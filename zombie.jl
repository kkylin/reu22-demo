# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Julia 1.7.1
#     language: julia
#     name: julia-1.7
# ---

# Modeling a Zombie Apocalypse
# ============================
#
# This is based on [a demonstration of how to solve a system of first order ODEs
# using SciPy](https://scipy-cookbook.readthedocs.io/items/Zombie_Apocalypse_ODEINT.html). Note that a Nth order equation can also be solved using
# SciPy by transforming it into [a system of first order
# equations](http://en.wikipedia.org/wiki/Ordinary_differential_equation#Reduction_to_a_first_order_system).
# In a this lighthearted example, a system of ODEs can be used to model a
# "zombie invasion", using the equations specified in [Munz et al.
# 2009](http://mysite.science.uottawa.ca/rsmith43/Zombies.pdf).
#
# The system is given as:

# \begin{align}
# S'(t) &= P - B\cdot S(t)\cdot Z(t) - d\cdot S(t)\\
# Z'(t) &= B\cdot S(t)\cdot Z(t) + G\cdot R(t) - A\cdot S(t)\cdot Z(t)\\
# R'(t) &= d\cdot S(t) + A\cdot S(t)\cdot Z(t) - G\cdot R(t)\\
# \end{align}

# with the following notation:
#
# *  S(t): the number of susceptible victims at time t
# *  Z(t): the number of zombies at time t
# *  R(t): the number of people "killed" by time t
# *  P: the population birth rate
# *  d: the chance of a natural death
# *  B: the chance the "zombie disease" is transmitted (an alive person becomes a zombie)
# *  G: the chance a dead person is resurrected into a zombie
# *  A: the chance a zombie is totally destroyed

# This involves solving a system of first order ODEs given by: $y'(t) =
# f(y(t), t)$ where $y(t) = [S(t), Z(t), R(t)]$.
#
# To run this notebook, you will also need to download the files [ode.jl](https://drive.google.com/file/d/1ZvmFf2iR-8fCJ_5EnP6MDDKQcET_8VxZ/view?usp=sharing) and [zombie.jl](https://drive.google.com/file/d/1PtELx0EIB2AjtSrdrtrwsxU-BHuxE9cY/view?usp=sharing).  Be sure to save them in the same directory as this notebook.

using PyPlot

include("rk2.jl")

# +
# zombie apocalypse modeling
function ZombieVectorField(;
        P = 0,       # birth rate
        d = 0.0001,  # natural death percent (per day)
        B = 0.0095,  # transmission percent  (per day)
        G = 0.0001,  # resurect percent (per day)
        A = 0.0001  # destroy percent  (per day)
        )

    # solve the system dy/dt = f(y, t)
    function(y, t)
        Si = y[1]
        Zi = y[2]
        Ri = y[3]
        # the model equations (see Munz et al. 2009)
        f0 = P - B*Si*Zi - d*Si
        f1 = B*Si*Zi + G*Ri - A*Si*Zi
        f2 = d*Si + A*Si*Zi - G*Ri
        return [f0, f1, f2]
    end
end

# initial conditions
S0 = 500.     	  # initial population
Z0 = 0        	  # initial zombie population
R0 = 0        	  # initial death population
y0 = [S0, Z0, R0] # initial condition vector
t  = 0:(5/1000):10 # time grid
# -

# solve the DEs
soln1 = rk2(ZombieVectorField(), y0, t)
S1 = soln1[:, 1]
Z1 = soln1[:, 2]
R1 = soln1[:, 3]

# solve the DEs
soln2 = rk2(ZombieVectorField(P=30), y0, t)
S2 = soln2[:, 1]
Z2 = soln2[:, 2]
R2 = soln2[:, 3]

plot(t,S1;label="Living")
plot(t,Z1;label="Zombie")
plot(t,R1;label="Dead")
plot(t,S2;label="Living'")
plot(t,Z2;label="Zombie'")
plot(t,R2;label="Dead'")
legend()
title("Zombie apocalypse")
xlabel("Days since outbreak")
ylabel("Population")


