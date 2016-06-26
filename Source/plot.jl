include("model.jl")
using PyPlot


#plot the states
function plotStates( Ref,     #reference trajectory
                     Xs       #states matrix
                   )
    Xs = Xs';
    Ref = Ref[1: nT-N, :];
    T = Ref[1: nT-N, 1];

    for i = 1: ns
        figure(i)
        plot(T, Ref[:, i+1], color="blue", linewidth=5.0);
        plot(T, Xs[:, i], color="red", linewidth=5.0, linestyle="--");
    end
end
