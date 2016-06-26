###########################################
#  This mudule is for 4/5Runge Kutta Method 
#  Author: Marvin
#  Date: 2015.12.20
#  Version: 1.2
#  Warning: init t must be 0!!


type OdeCoef
a::Array
b::Array
bs::Array
OdeCoef(a,b,bs)=new(a,b,bs)
end

const RunKuDP45 = OdeCoef(
		[0 0 0 0 0 0 0
		1/5 0 0 0 0 0 0
		3/40 9/40 0 0 0 0 0
		44/45 -56/15 32/9 0 0 0 0
		19372/6561 -25360/2187 64448/6561 -212/729 0 0 0
		9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0
		35/384 0 500/1113 125/192 -2187/6784 11/84 0],
		[35/384 0 500/1113 125/192 -2187/6784 11/84 0],
		[5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40],
		)

###########################################
#  This function is for calculate k 
#  k is a n by 7 matrix
#  Date: 2015.10.30
#  Version: 1.0
function Slope(f::Function,N::Integer,h::Number,y::Array,der::Bool)
k=Array(Number,N,7);
k=zeros(k)
cof=Array(Number,N,7);
cof=zeros(cof);
for i=1:7
	if der
		cof[:,i]=y+k*(RunKuDP45.a[i,:])'*h
	end
	k[:,i]=f(y+k*(RunKuDP45.a[i,:])'*h)
end

if der
	return [k cof]
else 
	return k
end
end

###########################################
#  This function is for calculate norm e 
#  e is a vector with n elements
#  Date: 2015.10.30
#  Version: 1.0
function StepSize(yn::Array,zn::Array,N::Integer,TOL::Number,ATOL::Number,y::Array)
e=Array(Number,N);
e=zeros(e);
for i=1:N
	e[i]=abs(yn[i]-zn[i])/(TOL*y[i]+ATOL)
end
norme=norm(e)
return norme
end

function ant_vec(vect::Array,N::Integer)
antvec_mat=Array(Float64,N,N)
k=zeros(antvec_mat)
for i=1:N
	antvec_mat[:,i]=vect[(i-1)*N+1:i*N]
end
return antvec_mat 
end

###########################################
#  This function is for calculate the Jacobian
#  return y which is a vector
#  Date: 2015.12.22
#  Version: 1.2
function Jacob(f::Function,jac_f::Function,N::Integer,h::Number,cof::Array)
jac_k=Array(Float64,N*N,7);
jac_k=zeros(jac_k);
for i=1:7
	jac_k[:,i]=vec(jac_f(cof[:,i])*(eye(N)+ant_vec(jac_k*(RunKuDP45.a[i,:])'*h,N)))
end
dy_2 = jac_k*(RunKuDP45.b)'*h
return ant_vec(dy_2,N) +eye(N)
end



###########################################
#  This function is for calculate the T current y(T)
#  return y which is a vector
#  Date: 2015.12.20
#  Version: 1.2
function RK45(f::Function,x0::Array,delT::Number,start_time::Number, hinit::Number,TOL::Number, ATOL::Number,jac_f::Function,der::Bool)
n = length(x0);
t=start_time
T=delT+start_time
#H=Array(Number,0);
h	=    hinit*sign(T); 
y	=    x0;
yn	=    y;
dy  	=    zeros(n,n);
dy_1	=    eye(n)
dyn 	=    eye(n)
dy_n_1  =    eye(n) 
en=Array(Number,0)
while 1==1
	t=t+h
	if abs(t)>= abs(T)
		h=T-t+h
	end
	while 1==1
		k=Slope(f,n,h,y,der);
		if der
			cof=k[:,8:14] 
			k=k[:,1:7]	
			dy_n_1 = Jacob(f,jac_f,n,h,cof) * dy_1
		end
		yn=k*(RunKuDP45.b)'*h+y;
		zn=k*(RunKuDP45.bs)'*h+y;
		norme=StepSize(yn,zn,n,TOL,ATOL,y)
		#H=[H;h]   Test step
		if norme <= 1
			break;
		end
        h=0.9*h/norme^(1/5)
#		em=yn-zn;
#		en=[en;em]
	end
	y = yn;
	if der 
		dy_1 = dy_n_1
	end
	if abs(t)>=abs(T)
		break;
	end		
	#Y=[Y;y]
end
# return en; Use this to determine that error is always 0.(for Question 2)
# return [y;H]; Test step
if der
	return dy_1
else
	return y;
end
end



