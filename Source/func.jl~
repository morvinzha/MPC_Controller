include("model.jl")

function initF( y,     # initial states
                x,     # states
                R      # periodical reference trajectory
              )
#    n = size(y, 1);
    x = reshape(x, n, N);
    x = x[1: ns, :];
    y = y[1: ns, :];
    F = vec(hcat(y, x) - R');
    return F;
end

function dF( y,     # initial states
             x,     # states
             d,     # derivation direction
             R      # periodical reference trajectory
           )
    h = 1.0e-6;
    if d == 'x'
        nx = size(x, 1);
        dFmtx = Array{Float64}(ns*(N+1), nx);
        for j = 1: nx
            t = x[:];
#            t[j] = t[j] + 10^(-100.0)*(-1+0im)^0.5;
            t[j] = t[j] + h;
#            dFmtx[:, j] = (initF(y, t, R))/10^(-100.0);
            dFmtx[:, j] = (initF(y, t, R) - initF(y, x, R))/h;
        end
    elseif d == 'y'
        ny = size(y, 1);
        dFmtx = Array{Float64}(ns*(N+1), ny);
        for j = 1: ny
            t = y[:];
            t[j] = t[j] + h;
            dFmtx[:, j] = (initF(t, x, R) - initF(y, x, R))/h;
        end
    end
    return dFmtx;
end

function initG( y,     # initial states
                x,     # states
                u,     # controls
              )
    n = size(y, 1);
    N = round(Int, size(x, 1)/n);
    m = round(Int, size(u, 1)/N);

    x = reshape(x, n, N);
    u = reshape(u, m, N);

    G = Array{Float64}(n, N);
    G[:, 1] = x[:, 1] - E(y, u[:, 1]);
    for i = 1: N-1
        G[:, i+1] = x[:, i+1] - E(x[:, i], u[:, i+1]);
    end
    G = vec(G);
    return G;
end

function E( x,     # states
            u,     # controls
          )
    h  = (Ref[nT, 1] - Ref[1, 1])/(nT - 1);
    k1 = myode(x, u);
    k2 = myode(x + h/2*k1, u);
    k3 = myode(x + h/2*k2, u);
    k4 = myode(x + h*k3, u);
    x = x + h*(1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4);
end

function dG( y,     # initial states
             x,     # states
             u,     # controls
             d,     # derivation direction
           )
    h = 1.0e-6;
    if d == 'x'
        nx = size(x, 1);
        dGmtx = Array{Float64}(n*N, nx);
        for j = 1: nx
            t = x[:];
            t[j] = t[j] + h;
            dGmtx[:, j] = (initG(y, t, u) - initG(y, x, u))/h;
        end
    elseif d == 'y'
        ny = size(y, 1);
        dGmtx = Array{Float64}(n*N, ny);
        for j = 1: ny
            t = y[:];
            t[j] = t[j] + h;
            dGmtx[:, j] = (initG(t, x, u) - initG(y, x, u))/h;
        end
    elseif d == 'u'
        nu = size(u, 1);
        dGmtx = Array{Float64}(n*N, nu);
        for j = 1: nu
            t = u[:];
            t[j] = t[j] + h;
            dGmtx[:, j] = (initG(y, x, t) - initG(y, x, u))/h;
        end
    end
    return dGmtx;
end
