module dvr

export exp_dvr
#############################################################################
function exp_dvr(N)

#Write DVR grid and position operators on grid#
phi = zeros((2*N+1))
X = zeros((2*N+1,2*N+1))
Y = zeros((2*N+1,2*N+1))
V6 = zeros((2*N+1,2*N+1))

for ii=1:(2*N+1)
	phi[ii] = ii*2.0*pi/(2*N+1)
	X[ii,ii] = cos(phi[ii])	
	Y[ii,ii] = sin(phi[ii])	
	V6[ii,ii] = cos(6.0*phi[ii])	
end

#Write kinetic energy matrix in DVR basis#
T = zeros((2*N+1,2*N+1))

for ii=1:(2*N+1)
	for jj=ii:(2*N+1)
		dind = ii-jj
		T[ii,jj] = cos(pi*dind/(2*N+1))/(2.0*sin(pi*dind/(2*N+1))^2)*(-1.0)^dind
		T[jj,ii] = T[ii,jj]
	end
	T[ii,ii] = N*(N+1)/3
end

#Write Lz operator in DVR basis#
Lz = zeros(ComplexF64,(2*N+1,2*N+1))

dphi = 2*pi/(2*N+1)

for ii=1:(2*N+1)
	for jj=1:(2*N+1)
		dum=0.0+0im
		for n=-N:N
			dum+=n*exp(im*n*(phi[ii]-phi[jj]))
		end
		Lz[ii,jj] = dum*dphi/(2*pi)
	end
end

return T,X,Y,Lz,V6
end
#############################################################################
end
