include("func.jl")
#include("plot.jl")
include("qpOASES.jl")
#using Convex

# optimal controls
ctrl = zeros(m, nT-N+1);

# initial state
x0 = zeros(n, 1);

# initial guess of states
y = zeros(n, 1);
x = zeros(n*N, 1);
u = zeros(m*N, 1);

X = Array{Float64}(n, 1);

# loop
for i = 1: nT-N

    #------ preparation with condensing -------#

    # shift the states and controls
    u = [u[m+1: m*N]; u[m*(N-1)+1: m*N] ];

    #x = [x[n+1: n*N]; x[n*(N-1)+1: n*N]];

    # precompute F, Fx, and Fy
    F = initF(y, x, Ref[i: i+N, 2: ns+1]);
    Fx = dF(y, x, 'x', Ref[i: i+N, 2: ns+1]);
    Fy = dF(y, x, 'y', Ref[i: i+N, 2: ns+1]);

    # precompute G, Gx, Gy and Gu
    G = initG(y, x, u);
    Gx = dG(y, x, u, 'x');
    Gy = dG(y, x, u, 'y');
    Gu = dG(y, x, u, 'u');

    # precompute Hessian matrix H, gradient matrix g, and constraints matrix A
    Ry = Fy - Fx/Gx*Gy;
    Ru = -Fx/Gx*Gu;
    R  = F - Fx/Gx*G;
    H  = [Ry Ru];
    g  = 2*R'*H;
    H  = 2*H'*H;
    A  = [eye(n) zeros(n, m*N)];

    #GradF = [Fx Fy zeros(ns*(N+1), m*N)];
    #H = 2 * GradF' * GradF;
    #g = 2*(F' * GradF)';
    #I = [zeros(n, n*N) eye(n) zeros(n, m*N)];
    #A = [Gx Gy Gu];
    #A = [A; I];

    #---------- receive measurement -----------#

    if i == 1
        X = x0;
    else
        x0 = E(x0, ctrl[:, i-1]);
        X = [X x0];
    end

    #------------ instant feedback ------------#

    # compute boundary vector bA with the receive measurement x0
    bA = x0 - y;

    #bA = [-G; x0-y];

    # solve SQP for one iteration to get a instant feedback
    nV = size(H, 1);
    nC = size(A, 1);
    delv = Array{Float64}(nV, 1);
    ##fval =
    QPOASES_JULIA(H, g, A, [], [], bA, bA, delv, Int64(nV), Int64(nC));

    # generate iterates of the form
    y = y + delv[1: n];
    u = u + delv[n+1: n+m*N];

    #x = x + delv[1:n*N];
    #y = y + delv[n*N+1: n*N+n];
    #u = u + delv[n*N+n+1: n*N+n+m*N];

    # send first control to the real process
    ctrl[:, i] = u[1: m, 1];

end

# plot states trajectories
#plotStates( Ref, X );
