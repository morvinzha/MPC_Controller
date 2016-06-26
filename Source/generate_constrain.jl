include("ODE.JL")
include("dynamics.jl")


function NDforward(f::Function,x0)
h=(eps())^(2/3)
n=length(x0)
E=eye(n)
df=zeros(n,n)
for i=1:n
	df[:,i]=(f(x0+h*E[:,i])-f(x0-h*E[:,i]))/(2*h)
end
return df
end


################################
#   	Multiple shooting
#
#
################################

################################
#	RK45 integrator
#	Input				--				Output
#	y									yn
#	u
#	start_time
#	sample_time(global)
#	RK45_init,Tol,ATol(global)
#	dyn()(dynamics.jl)
#	jac_dyn(dynamics.jl)
function next_state(y,u,start_time)
#	current_state = RK45(dyn,[y;u],sample_time,start_time,RK45_hinit,RK45_Tol,RK45_ATol,jac_dyn,false) #
#	current_state = RK45(dyn_uav_f,[y;u],sample_time,start_time,RK45_hinit,RK45_Tol,RK45_ATol,jac_dyn_f,false)
	current_state = RK45(dyn_f,[y;u],sample_time,start_time,RK45_hinit,RK45_Tol,RK45_ATol,jac_dyn_f,false)
	yn	  = current_state[1:dim_sample_x]
	return yn
end

################################
#	Derivate of RK45 integrator to its init value
#	Input			    --				Output
#	y									dRK45/dy
#	u
#	start_time
#	sample_time(global)
#	RK45_init,Tol,ATol(global)
#	dyn()(dynamics.jl)
#	jac_dyn(dynamics.jl)
function sample_jac_xu(y,u,start_time)
#	jac = RK45(dyn,[y;u],sample_time,start_time,RK45_hinit,RK45_Tol,RK45_ATol,jac_dyn,true) #
	jac = RK45(dyn_f,[y;u],sample_time,start_time,RK45_hinit,RK45_Tol,RK45_ATol,jac_dyn_f,true)
#function G(x)
#	return RK45(dyn_uav_f,x,sample_time,start_time,RK45_hinit,RK45_Tol,RK45_ATol,jac_dyn,false)
#end
#	jac = NDforward(G,[y;u])
	return jac
end


################################
#	Build the first term of A_dyn by compose each node dRK45/dy
#	Input			    --				Output
#	A_x:[dRK45/dy dRK45/dy dRK45/dy]	A_x: The first term of A_dyn 
#
function column_exchange(A_x::Array)
	A_sx  = Array(Float64,dim_sample_x*times_of_shooting,0)
	A_su  = Array(Float64,dim_sample_x*times_of_shooting,0)
	A_x  = A_x[1:dim_sample_x,:]						# delete control part
	index_sample_x = 1
	index_sample_u = dim_sample_x + 1  							
	for i=1:times_of_shooting
		A_xmat = A_x[:,index_sample_x:index_sample_x+dim_sample_x-1]
		A_xmat = [zeros((i-1)*dim_sample_x,dim_sample_x);A_xmat;zeros(nC_dyn-i*dim_sample_x,dim_sample_x)]
		A_sx=[A_sx A_xmat]
		A_umat = A_x[:,index_sample_u:index_sample_u+dim_sample_u-1]
		A_umat = [zeros((i-1)*dim_sample_x,dim_sample_u);A_umat;zeros(nC_dyn-i*dim_sample_x,dim_sample_u)]
		A_su   = [A_su A_umat]
		index_sample_x = index_sample_x+dim_sample_xu
		index_sample_u = index_sample_u+dim_sample_xu
	end
	A_s=[A_sx zeros(dim_sample_x*times_of_shooting,dim_sample_x) A_su]
#	A_s0 = [A_x[:,1:dim_sample_x];zeros(dim_sample_x,dim_sample_x)]
#	A_s1 = [zeros(dim_sample_x,dim_sample_x);A_x[:,dim_sample_xu+1:dim_sample_xu+dim_sample_x]]
#	A_s2 = zeros(nC_dyn,dim_sample_x)
#	A_c1 = [A_x[:,dim_sample_x+1:dim_sample_xu];zeros(dim_sample_x,dim_sample_u)] 
#	A_c2 = [zeros(dim_sample_x,dim_sample_u); A_x[:,dim_sample_xu+dim_sample_x+1:2*dim_sample_xu]]
#	A_x  = [A_s0 A_s1 A_s2 A_c1 A_c2]
	return A_s
end

################################
#	C(S) in G(S) = C(S) - S = 0
#	Input			    --				Output
#	
#
#
#
function MultiShoX_constrain()				#MultiShox(xOpt)=0
	state_mat = Array(Float64,dim_sample_x,times_of_shooting)
	state_mat = zeros(state_mat)
	index_state		= 1	
	index_control	= nV-number_control+1
	for i=1:times_of_shooting
		state_mat[:,i] = next_state(S[index_state:index_state+dim_sample_x-1],S[index_control:index_control+dim_sample_u-1],sample_time * (i-1+times_of_shift))
#		state_mat[:,i] = next_state(S[index_state:index_state+dim_sample_x-1],S[index_control:index_control+dim_sample_u-1],sample_time* (i-1))
		index_state    = index_state + dim_sample_x	
		index_control  = index_control + dim_sample_u
	end	
	state_varible  = vec(state_mat)
	return state_varible
end

################################
#	dC(s) in dG(S) = dC(S) - dS = 0
#	Input			    --				Output
#
#
#
#
function MultiShoX_dG()
	state_mat = Array(Number,nC_dyn+2,nV)
	state_mat = zeros(state_mat)
	index_state = 1
	index_control = dim_sample_x * (times_of_shooting + 1)+1
	state_mat = sample_jac_xu(S[index_state:index_state+dim_sample_x-1],S[index_control:index_control+dim_sample_u-1],times_of_shift*sample_time)
#	state_mat = sample_jac_xu(S[index_state:index_state+dim_sample_x-1],S[index_control:index_control+dim_sample_u-1],0)
	for i=1:times_of_shooting - 1
		index_state   = index_state + dim_sample_x
		index_control = index_control + dim_sample_u
#		state_mat = [state_mat sample_jac_xu(S[index_state:index_state+dim_sample_x-1],S[index_control:index_control+dim_sample_u-1], sample_time*(i+times_of_shift))]
		state_mat = [state_mat sample_jac_xu(S[index_state:index_state+dim_sample_x-1],S[index_control:index_control+dim_sample_u-1], sample_time*i)]
	end
#	state_mat = [state_mat sample_jac_xu(state_1,control_2,0.5)]
	dG 		  = column_exchange(state_mat)
	return dG
end

################################
#   	generate constrain
#
#
################################


################################
#
#
function generate_A_dyn()
	dG1_s	= MultiShoX_dG()
	dG2_s	= [zeros(nC_dyn,dim_sample_x) eye(nC_dyn) zeros(nC_dyn,times_of_shooting * dim_sample_u)]
	A_dyn	= dG1_s - dG2_s
	return A_dyn
end


################################
#
#
function generate_lubA_dyn()
	state_vec=Array(Float64,0)
	for i=1:times_of_shooting	
		state_vec=[state_vec;S[i*dim_sample_x+1:(i+1)*dim_sample_x]]
	end
	lbA_dyn=-MultiShoX_constrain()+state_vec
#	lbA_dyn=-MultiShoX_constrain()+[state_1;state_2]
	return lbA_dyn
end

################################
#
#
function generate_A_eq()
	A_eq=zeros(nC_eq,nV)
	A_eq[1:position_freedom+speed_freedom,1:position_freedom+speed_freedom]=eye(position_freedom+speed_freedom)	
	A_eq[position_freedom+speed_freedom+1:nC_eq,nV-number_control-dim_sample_x+1:nV-number_control-1]=eye(position_freedom+speed_freedom)
	return A_eq
end

################################
#
#
function generate_lubA_eq(A_eq)
	b_eq=zeros(nC_eq)
	b_eq=[start_lbound;end_ubound]
	lbA_eq=b_eq-A_eq*S
	return lbA_eq
end

################################
#
#
function generate_A_neq()
#	A_neq=zeros(nC_neq,nV)
#	A_neq[1:number_control,nV-number_control+1:nV]=eye(number_control)
	A_neq=eye(nV)
	return A_neq
end

################################
#
#
function generate_lbA_neq(A_neq)
	lb_neq=zeros(nC_neq)
	lb_neq=[state_lbound;control_lbound]
	lbA_neq=-A_neq*S+lb_neq
	return lbA_neq
end

################################
#
#
function generate_ubA_neq(A_neq)
	ub_neq=zeros(nC_neq)
	ub_neq=[state_ubound;control_ubound]
	ubA_neq=ub_neq-A_neq*S
	return ubA_neq
end

function fobs(s,obs_sample)
	n=length(s)
	norm_2=norm(s-obs_sample)
	return norm_2^2
end

function dfobs(s,obs_sample)
	n=length(s)
	norm_2=Array(Float64,0)
	r=s-obs_sample;
	for i=1:n
		norm_2=[norm_2;2*abs(obs_sample[i]-s[i])]
	end
	return norm_2
end

function MultiShotObs_fobs(S,obs_sample)
	ubA_obs=Array(Float64,0)
	for i=0:times_of_shooting
		ubA_obs=[ubA_obs;fobs(S[i*dim_sample_x+1:i*dim_sample_x+dim_sample_x_r],obs_sample)]
	end
	return ubA_obs
end

function MultiShotObs_dfobs(S,obs_sample)
	A_obs_row=Array(Float64,0)
	A_obs=Array(Float64,0,length(S)-times_of_shooting*dim_sample_u)	
	for i=0:times_of_shooting
		A_obs_row=[dfobs(S[i*dim_sample_x+1:i*dim_sample_x+dim_sample_x_r],obs_sample);zeros(dim_sample_x-dim_sample_x_r,1)]
		A_obs=[A_obs;[zeros(1,i*dim_sample_x) A_obs_row' zeros(1,(times_of_shooting-i)*dim_sample_x)]]
	end
	A_obs=[A_obs zeros(times_of_shooting+1,times_of_shooting*dim_sample_u)]
	return A_obs
end


function generate_A_obs(S,obs_sample)
	A_obs=MultiShotObs_dfobs(S,obs_sample)
	return A_obs
end


function generate_ubA_obs(S,obs_sample)
	ubA_obs=MultiShotObs_fobs(S,obs_sample)
	ubA_neq=fill(obs_ub,times_of_shooting+1)-ubA_obs
	return ubA_neq
end

function generate_lbA_obs(S,obs_sample)
	lbA_obs=MultiShotObs_fobs(S,obs_sample)
	lbA_neq=-lbA_obs+fill(obs_lb,times_of_shooting+1)
	return lbA_neq
end
