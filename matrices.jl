module matrices

using LinearAlgebra

export kinetic,Xoperator,Yoperator,lz_operator,V6_operator
#################################################################
function kinetic(mmax)

matrix = zeros((2*mmax+1),(2*mmax+1))

k=0
for m=-mmax:mmax
	k+=1
	matrix[k,k] = m*m
end	

return matrix
end
#################################################################
function Xoperator(mmax)

diagonal = zeros((2*mmax+1))
subdiagonal = 0.5*ones((2*mmax))

matrix=Tridiagonal(subdiagonal,diagonal,subdiagonal)

matrix_out = zeros(Float64,(2*mmax+1,2*mmax+1))
for i1=1:2mmax+1
for i2=1:2mmax+1
	matrix_out[i1,i2]= matrix[i1,i2]
end
end

return matrix_out
end
#################################################################
function Yoperator(mmax)

diagonal = zeros((2*mmax+1))
upper_subdiagonal = 0.5*ones((2*mmax))
lower_subdiagonal = -0.5*ones((2*mmax))

matrix=Tridiagonal(lower_subdiagonal,diagonal,upper_subdiagonal)

matrix_out = zeros(Float64,(2*mmax+1,2*mmax+1))
for i1=1:2mmax+1
for i2=1:2mmax+1
	matrix_out[i1,i2]= matrix[i1,i2]
end
end

return matrix_out
end
#################################################################
function lz_operator(mmax)

matrix = zeros((2*mmax+1),(2*mmax+1))

k=0
for m=-mmax:mmax
	k+=1
	matrix[k,k] = m
end	

return matrix
end


#################################################################
function V6_operator(mmax)

matrix_out = zeros(Float64,(2*mmax+1,2*mmax+1))
for i1=1:2mmax+1
for i2=1:2mmax+1
	if abs(i1-i2)==6
		matrix_out[i1,i2]= 0.5
	end
end
end

return matrix_out
end

#################################################################
end 
