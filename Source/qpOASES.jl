#----------- QP Solver Interface -----------#
#----------- min  x'Hx + g'x ---------------#
#----------- s.t. lb  <=  x <= ub ----------#
#---------------- lbA <= Ax <= ubA ---------#
#------------ Author Xuhui Feng ------------#
#------------- Date 11/20/2015 -------------#

function QPOASES_JULIA( H::Array{Float64},         # Hessian matrix
                        g::Array{Float64},         # Gradient vector
                        A::Array{Float64},         # Constraint matrix
                        lb::Array,                 # Lower bound vector (on variables)
                        ub::Array,                 # Upper bound vector (on variables)
                        lbA::Array,                # Lower constraints' bound vector
                        ubA::Array,                # Upper constraints' bound vector
                        xOpt::Array{Float64},      # Optimal solution
                        nV::Int64,                 # Number of variables
                        nC::Int64                  # Number of constraints
                      )
    H = vec(H);
    A = vec(A');

    for i=1:length(lb)
	if ub[i]==Inf && ub[i] ==NaN
		ub[i] = C_NULL;
	elseif lb[i]==Inf && lb[i] ==NaN
		lb[i] = C_NULL;
	end
    end

    for i=1:length(lbA)
	if ubA[i]==Inf && ubA[i] ==NaN
		ubA[i] = C_NULL;
	elseif lbA[i]==Inf && lbA[i] ==NaN
		lbA[i] = C_NULL;
	end
    end


    if lb == []
        lb = C_NULL;
    end

    if ub == []
        ub = C_NULL;
    end

    if lbA == []
        lbA = C_NULL;
    end

    if ubA == []
        ubA = C_NULL;
    end



    fval = ccall((:qpOASES, "/home/marvin/Julia_linux/nlMPC_Realtime/libqpOASES.so"),
                 (Float64),
                 (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Int64),
                 H, g, A, lb, ub, lbA, ubA, xOpt, nV, nC);

    return fval;
end
