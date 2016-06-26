using Gadfly

function plot_all()
	plot(x=result[:,1],y=result[:,2],Geom.point,Geom.line)
	plot(x=result[:,1],y=result[:,3],Geom.point,Geom.line)
	plot(x=result[:,1],y=result[:,4],Geom.point,Geom.line)
	plot(x=result[:,1],y=result[:,5],Geom.point,Geom.line)
end

function writedata()
f=open("Data1.dat",false,true,true,false,true)
writedlm(f,result[:,1:2])
flush(f)
close(f)

f=open("Data2.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,3]])
flush(f)
close(f)

f=open("Data3.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,4]])
flush(f)
close(f)

f=open("Data4.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,5]])
flush(f)
close(f)

f=open("Data5.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,6]])
flush(f)
close(f)

f=open("Data6.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,8]])
flush(f)
close(f)

f=open("Data7.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,9]])
flush(f)
close(f)

f=open("Data8.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,3]])
flush(f)
close(f)

f=open("DataXY.dat",false,true,true,false,true)
writedlm(f,[result[:,2] result[:,3]])
flush(f)
close(f)

f=open("DataC.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,14]])
flush(f)
close(f)
f=open("DataC.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,15]])
flush(f)
close(f)
f=open("DataC.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,16]])
flush(f)
close(f)
f=open("DataC.dat",false,true,true,false,true)
writedlm(f,[result[:,1] result[:,17]])
flush(f)
close(f)

end


