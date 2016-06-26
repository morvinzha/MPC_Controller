include("qpOASES.jl")
include("generate_constrain.jl")
include("dynamics.jl")
include("MPCplot.jl")

################################
#    Init Globle Constant
#
#
################################
# constant for MPC controller
times_of_shooting = 10									# multiple shooting
sample_time       = 1/times_of_shooting					# discreatimilize Real time = sample_t*T

# constant for Dynamic model								# Rotate inercial
#const position_freedom		= 1;
const position_freedom	= 3;
#const position_freedom  = 6;

#const speed_freedom		= 1;
const speed_freedom		= 2;
#const speed_freedom		= 6;

#const dim_sample_x 		= 3;									# dim of sample x;
const dim_sample_x 		= 6;#
#const dim_sample_x 		= 13;

#const dim_sample_u			= 1;
const dim_sample_u		= 2;
#const dim_sample_u		= 4;

const dim_sample_xu		= dim_sample_x + dim_sample_u

const dim_sample_x_r	= 2;	# For generate obstacles
#const dim_sample_x_r	= 3;

global number_control	= dim_sample_u * times_of_shooting

#const T_position			= 3
const T_position		= 6
#const T_position		= 13

# constant only for RK45
const RK45_Tol   = 0.0001;
const RK45_ATol  = 0.000001;
const RK45_hinit = 0.001;

# constant for qpOASES Sover
global nV     = (times_of_shooting + 1) * dim_sample_xu - dim_sample_u  	 	# Number of variables  (times_of_shooting + 1)* dim of sample [x;u] - end u
global nC_dyn = times_of_shooting * dim_sample_x         	    				# Number of constraints  times_of_shooting * dim of sample x;
global nC_neq = nV #times_of_shooting * dim_sample_u								# nonlinear ineq contrains   dim of u +
const nC_eq   = 2 * (position_freedom + speed_freedom) 							# linear eq contrains	 start + end
global nC_obs  = times_of_shooting+1
## constant for obj function "minize x'Hx+g'x"
H							= zeros(nV,nV); 	        						# Hessian matrix
#H[T_position,T_position]	= 1													
g							= zeros(nV);             							# Gradient vector
g[T_position]=1.0;

## constant for constrains
### constant for start state
#const r_0 = 0.0
#const v_0 = 0.0
 rx_0    = 0.0
 ry_0    = 0.0
 Phi_0   = pi/4
 V_0     = 20.0
 Omega_0 = 0.0
#rx_0	=	0;
#ry_0	= 	0;
#rz_0	= 	0;
#vx_0 = 	0;
#vy_0 = 	0;
#vz_0 = 	0;
#roll_0	=	0;
#pitch_0=	0;
#yaw_0	=	0;
#droll_0=	0;
#dpitch_0=	0;
#dyaw_0=	0;
### constant for state
 rx_min    = 0
 rx_max    = 80

 ry_min    = 0
 ry_max    = 80

 V_min   = -28.0
 V_max   =  28.0

 Phi_min = -pi
 Phi_max = pi

 Omega_min = -pi/4
 Omega_max = pi/4

T_min = 0
T_max = 30

#rx_min	=	0;
#ry_min	= 	0;
#rz_min	= 	0;

#rx_max	=	50;
#ry_max	= 	50;
#rz_max	= 	50;

#vx_min = 	-20;
#vy_min = 	-20;
#vz_min = 	-20;

#vx_max = 	20;
#vy_max = 	20;
#vz_max = 	20;

#roll_min	= -1*pi/3;;
#pitch_min	= -1*pi/3;
#yaw_min		= 0	

#roll_max	=1*pi/3	;
#pitch_max	=1*pi/3;
#yaw_max		= 2*pi;

#droll_min	=-1*pi/6; 
#dpitch_min	=-1*pi/6;
#dyaw_min	=-1*pi/6;

#droll_max	=1*pi/6;
#dpitch_max	=1*pi/6;
#dyaw_max	=1*pi/6;

### constant for end state
#const r_T = 1.50
#const v_T = 0.0
 rx_T    = 50
 ry_T    = 50
 V_T     = 10.0
 Phi_T	  = pi/4
 Omega_T = 0.0

#rx_T	=	0;
#ry_T	= 	0;
#rz_T	= 	20;
#vx_T = 	0;
#vy_T = 	0;
#vz_T = 	0;
#roll_T	=	0;
#pitch_T=	0;
#yaw_T	=	0;
#droll_T=	0;
#dpitch_T=	0;
#dyaw_T=	0;
### constant for control 
 a_min     = -3.0
 Ang_a_min = -0.1
 a_max     =  3.0
 Ang_a_max =  0.1
#w1_min		=  0
#w2_min		=  0
#w3_min		=  0
#w4_min		=  0
#w1_max		=  1000
#w2_max		=  1000
#w3_max		=  1000
#w4_max		=  1000

### constant for obstacles
obs_lb=4
obs_ub=6400
obs_x_position = 32
obs_y_position = 30
#obs_z_position = 30
global obs=[obs_x_position;obs_y_position]
#global obs=[obs_x_position;obs_y_position;obs_z_position]
obs_sample=obs
################################
#    Init Globle Peremeter
#
#
################################
global start_lbound	= [rx_0;ry_0;V_0;Phi_0;Omega_0]
#start_lbound	= [r_0;v_0]
#global start_lbound	= [rx_0;ry_0;rz_0;vx_0;vy_0;vz_0;roll_0;pitch_0;yaw_0;droll_0;dpitch_0;dyaw_0]

global end_ubound = [rx_T;ry_T;V_T;Phi_T;Omega_T]
#global end_ubound	= [rx_T;ry_T;rz_T;vx_T;vy_T;vz_T;roll_T;pitch_T;yaw_T;droll_T;dpitch_T;dyaw_T]
#global end_ubound	= [r_T;v_T]
global control_lbound=Array(Float64,0)
global control_ubound=Array(Float64,0)
for i=1:times_of_shooting
	control_lbound	= [control_lbound;a_min;Ang_a_min]
	control_ubound	= [control_ubound;a_max;Ang_a_max]
#	control_lbound	= [control_lbound;a_min]
#	control_ubound	= [control_ubound;a_max]
#	control_lbound	= [control_lbound;w1_min;w2_min;w3_min;w4_min]
#	control_ubound	= [control_ubound;w1_max;w2_max;w3_max;w4_max]
end

global state_lbound = Array(Float64,0) 
global state_ubound = Array(Float64,0)
for i=1:times_of_shooting+1
	state_lbound 	= [state_lbound ;rx_min;ry_min;V_min;Phi_min;Omega_min;T_min]
	state_ubound 	= [state_ubound ;rx_max;ry_max;V_max;Phi_max;Omega_max;T_max]
#	state_lbound 	= [state_lbound ;rx_min;ry_min;rz_min;vx_min;vy_min;vz_min;roll_min;pitch_min;yaw_min;droll_min;dpitch_min;dyaw_min;T_min]
#	state_ubound 	= [state_ubound ;rx_max;ry_max;rz_max;vx_max;vy_max;vz_max;roll_max;pitch_max;yaw_max;droll_max;dpitch_max;dyaw_max;T_max]
end


global xOpt		= zeros(nV)  									# Optimal solution	
global A		= zeros(nC_dyn+nC_eq+nC_neq,nV)

global lbA		= zeros(nC_dyn+nC_eq+nC_neq+nC_obs)               		# Lower constraints' bound vector
global ubA		= zeros(nC_dyn+nC_eq+nC_neq+nC_obs)					# Upper constraints' bound vector
                 											
global lb	             	   									# Lower bound vector (on variables)
global ub                										# Upper bound vector (on variables)		
global fval

################################
#    Init xOpt
#
#
################################														
xOpt[1:dim_sample_x-1] 				= start_lbound
for i=1:times_of_shooting+1
	xOpt[i*T_position] = 1
end
xOpt[nV-number_control+1:nV]		= control_ubound
#xOpt[times_of_shooting*+dim_sample_x+1:(times_of_shooting+1)*+dim_sample_x-1] 								= end_ubound

#function main()
global times_of_shift = 0;
result=Array(Float64,0,dim_sample_xu)
result=[start_lbound;0;zeros(dim_sample_u)]'# 0 is for parameter T
global S = zeros(nV)
 
for j = 1 : times_of_shooting-1
	lppoo = 2
	if times_of_shift < 4
		lppoo = 5
	end
	for i = 1 : lppoo
		S = S + xOpt
		A_dyn	= generate_A_dyn();						#generate contrains of dynamic system
		lbA_dyn	= generate_lubA_dyn();
		ubA_dyn	= lbA_dyn;
	
		A_eq	= generate_A_eq();						#generate linear equallity contrains of parameter T and start and end condition
										
		lbA_eq	= generate_lubA_eq(A_eq);				#generate b in Ads = (As - b) =
		ubA_eq	= lbA_eq;

		A_neq 	= generate_A_neq();						#generate linear inequlity contrains of s
		ubA_neq	= generate_ubA_neq(A_neq);				#generate ub in Ads < (As - ub)
		lbA_neq	= generate_lbA_neq(A_neq);				#generate lb in Ads > (As - lb)
	
		A_obs   = generate_A_obs(S,obs_sample);
		ubA_obs	= generate_ubA_obs(S,obs_sample);				
		lbA_obs	= generate_lbA_obs(S,obs_sample);		
		#conpose
		A=[A_dyn;A_eq;A_neq;A_obs];
		lbA=[lbA_dyn;lbA_eq;lbA_neq;lbA_obs];
		ubA=[ubA_dyn;ubA_eq;ubA_neq;ubA_obs];
		#A_single=[A;-A]
		#ubA_single=[ubA;-lbA]
	
		fval	= QPOASES_JULIA(H,g,A,[],[],lbA,ubA,xOpt,nV,nC_dyn+nC_eq+nC_neq)
#		fval	= QPOASES_JULIA(H,g,A_single,[],[],[],ubA_single,xOpt,nV,nC_dyn+nC_eq+nC_neq)
	end


times_of_shift 	+= 1
Accept_position = S[dim_sample_x+1:2*dim_sample_x]+xOpt[dim_sample_x+1:2*dim_sample_x]
Accept_control  = S[nV-number_control+1:nV-number_control+dim_sample_u] + xOpt[nV-number_control+1:nV-number_control+dim_sample_u]
result			= [result;[Accept_position' Accept_control']]


xOpt= [xOpt[1:(times_of_shooting+1)*dim_sample_x];xOpt[(times_of_shooting+1)*dim_sample_x+dim_sample_u+1:nV]]
S= [S[1:(times_of_shooting+1)*dim_sample_x];S[(times_of_shooting+1)*dim_sample_x+dim_sample_u+1:nV]]

xOpt= xOpt[dim_sample_x+1:length(xOpt)]
S   = S[dim_sample_x+1:length(S)]

times_of_shooting -= 1
number_control = number_control-dim_sample_u

nV	= nV - dim_sample_xu
nC_dyn	= nC_dyn - dim_sample_x
nC_neq = nC_neq - dim_sample_xu
nC_obs  = nC_obs - 1
start_lbound = Accept_position[1:dim_sample_x-1]
state_lbound = state_lbound[dim_sample_x+1:length(state_lbound)]
state_ubound = state_ubound[dim_sample_x+1:length(state_ubound)]
control_lbound = control_lbound[dim_sample_u+1:length(control_lbound)]
control_ubound = control_ubound[dim_sample_u+1:length(control_ubound)]

end
S = S + xOpt
result[1,T_position]=result[2,T_position] 																	# init time modified
result = [result;[(S[dim_sample_x+1:2*dim_sample_x])' (S[2*dim_sample_x+1:2*dim_sample_x+dim_sample_u])']]		# 
temp = result[:,dim_sample_x+1:dim_sample_x+dim_sample_u]
result[:,dim_sample_x+1:dim_sample_x+dim_sample_u] = [temp[2:size(temp,1),:];zeros(1,dim_sample_u)]

result = [result[2,T_position]*collect(0:sample_time:1) result]# time axis
Graph=plot_all()

plot(x=result[:,2],y=result[:,3],Geom.point,Geom.line)
#end
#ttmp=writedata()




