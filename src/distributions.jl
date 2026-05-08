module distributions

using LinearAlgebra

export distro_dvr,distro_complex
########################################################################################
function distro_dvr(N,phi) 

Nphi = length(phi)
#phi_dvr = [ii*2.0*pi/(2*N+1) for ii=1:(2*N+1)]
#phi = [ii*2.0*pi/Nphi for ii=1:Nphi] 

#D = zeros(Float64,(2*N+1,2*N+1,Nphi))
#U = zeros(ComplexF64,(Nphi,2*N+1))
#
#O = distro_complex(N,Nphi)
#
#
#for ip=1:Nphi
#	for i1=1:(2*N+1)
#		for i2=1:(2*N+1)
#			for m1=-N:N
#			for m2=-N:N
#				D[i1,i2,ip]+=real(O[N+m1+1,N+m2+1,ip]*exp(im*(m1*phi_dvr[i1]-m2*phi_dvr[i2]))/(2*pi))
#			end
#			end
#		end
#
##		dum=0.0+0im
##		for n=1:N-1
##			dum+=exp(im*n*(phi_dvr[i1]-phi[ip]))
##		end
##		U[ip,i1] = dum/(2*pi)
#	end
#end

#for ip=1:Nphi
#for i1=1:(2*N+1)
#for i2=1:(2*N+1)
#	D[i1,i2,ip] = real(conj(U[ip,i1])*U[ip,i2])
#end
#end
#end
#For grid at the DVR points#
D = zeros(Float64,(2*N+1,2*N+1,Nphi))
dphi = 2.0*pi/(2*N+1)
for ip = 1:Nphi 
	D[ip,ip,ip] = 1.0/dphi
end

return D
end
########################################################################################
function distro_complex(mmax,phi)

Nphi = length(phi)
D = zeros(ComplexF64,(2*mmax+1,2*mmax+1,Nphi))

for m1=-mmax:mmax
for m2=-mmax:mmax
	for ip=1:Nphi
		D[mmax+m1+1,mmax+m2+1,ip] = exp(im*(m2-m1)*phi[ip])/(2*pi)
	end
end
end

return D
end
########################################################################################
end
