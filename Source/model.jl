# reference trajectory
Ref = [0                  0                   0                   0                   ;
       0.250000000000000	0.00160595110888526	0.00808919130852264	1.51910671720685e-05;
       0.500000000000000	0.00634115024044085	0.0311587762664934	6.48988701534089e-05;
       0.750000000000000	0.0140417030283628	0.0668613955979837	0.000114971554713644;
       1                  0.0245432784797772	0.112850567454222	  0.000144355610175376;
       1.25000000000000	  0.0376807138590316	0.166780524268783	  0.000156951290793155;
       1.50000000000000	  0.0532861864229238	0.226306049357070	  0.000236459755131975;
       1.75000000000000	  0.0711963331970691	0.289082315202614	  0.000332781765898815;
       2	                0.0912574445300031	0.352764722946875	  0.000106221455316507;
       2.25000000000000	  0.113288001824743	  0.415008729605037	  5.13818292431731e-05;
       2.50000000000000	  0.137085314961041	  0.473469669517307	  0.00129797296116785 ;
       2.75000000000000	  0.162587748967510	  0.525802716262629	  0.000743941568910252;
       3	                0.189719817197608	  0.569662176211075	 -0.00429695427615495 ;
       3.25000000000000	  0.217786052957407	  0.602702983709079	  0.00208883717773625 ;
       3.50000000000000	  0.246608534096533	  0.622570870426288	  0.0202869090815436  ;
       3.75000000000000	  0.278412783847492	  0.626931189493662	 -0.0214309197713082  ;
       4	                0.311416463311287	  0.613344678643210	 -0.0745170676116542  ;
       4.25000000000000	  0.335912031511755	  0.579356186549952	  0.147809048920163   ;
       4.50000000000000	  0.364932986844674	  0.522432250544446	  0.245322794866297   ;
       4.75000000000000	  0.432918444115697	  0.431173876614494	 -0.834727338656262   ;
       5	                0.500000000000000	  0.300000000000000	 -1.57079632679490   ];

# number of sampling
nT = size(Ref, 1);

# MPC shift time horizon
N = 4;

# counts of states and controls
n  = 6;
ns = 3;
m  = 3;

# ordinary differential equations
function myode( x,     # states
                u,     # controls
              )
    M = 2.0;
    g = 9.81;
    gama_x = 0.1;
    gama_y = 0.1;
    gama_phi = 0.1;
    m = 1.0;
    l = 0.1;
    I = 0.02;

    u_x = u[1];
    u_y = u[2];
    u_phi = u[3];

    A = [     M + m              0          m * l * cos( x[3] );
                0              M + m        m * l * sin( x[3] );
           m * l * cos( x[3] )  m * l * sin( x[3] )    I      ];
    b = [   m * l * x[6]^2 * sin( x[3] ) + u_x + gama_x * x[4]                ;
          - m * l * x[6]^2 * cos( x[3] ) - ( M + m ) * g + u_y + gama_y * x[5];
          - m * g * l * sin( x[3] ) + u_phi - gama_phi * x[6]                ];

    dx = Array{Float64}(6, 1);
    dx[1] = x[4];
    dx[2] = x[5];
    dx[3] = x[6];
    dx[4:6] = A\b;

    return dx;
end
