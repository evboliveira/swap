module expectations
  
using ITensors,ITensorMPS
using LinearAlgebra


export vN_entropy,polarization,correlation,get_pphi,ang_correlation,binder,binder2
################################################################################
function vN_entropy(wf,mbond)

if length(wf) == 2
        orthogonalize!(wf, 2)
        U,S,V = svd(wf[2], (siteind(wf,2)))
else
        orthogonalize!(wf, mbond)
        U,S,V = svd(wf[mbond], (linkind(wf, mbond-1), siteind(wf,mbond)))
end
SvN = 0.0
renyi = 0.0
schmidtvalues=zeros(dim(S, 1))
for n=1:dim(S, 1)
        p = S[n,n]^2
        schmidtvalues[n]=p
        SvN -= p * log(p)
	renyi += p^2
end
renyi=-log(renyi)

return SvN,renyi,schmidtvalues
end
################################################################################
function polarization(wf,Nsites,Nbasis,evod,Xmat,Ymat)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis

include("operators.jl")

mux = expect(wf,"X")
muy = expect(wf,"Y")

if evod == "all" 
	muy = expect(wf,"Ycomp")
end
	
return real(sum(mux)),real(sum(muy))
end
################################################################################
function correlation(wf,Nsites,Nbasis,evod,Xmat,Ymat)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis

include("operators.jl")

dumX = correlation_matrix(wf,"X","X")
dumY = correlation_matrix(wf,"Y","Y")
if evod == "all"
	dumY = -dumY
end

Xcorr=0.0
Ycorr=0.0
for ii=1:Nsites-1
	Xcorr+=dumX[ii,ii+1]
	Ycorr+=dumY[ii,ii+1]
end

return real(Xcorr),real(Ycorr)
end
################################################################################
function ang_correlation(wf,Nsites,Nbasis,evod,Lzmat)

global Lz = Lzmat
global Nspec = Nbasis

include("operators.jl")

dumLz = correlation_matrix(wf,"Lz","Lz")

if evod == "all_real" || evod == "dvr"
	dumLz = -dumLz
end

Lzcorr=0.0
for ii=1:Nsites-1
	Lzcorr+=dumLz[ii,ii+1]
end

return real(Lzcorr)
end
################################################################################
function get_pphi(psi,Nsites,Nphi,Nbasis,Dphi)

global D = Dphi
global Nspec = Nbasis

phi_dist = zeros(Float64,(Nsites,Nphi))
for ip=1:Nphi
	dum = expect(psi,string("D",ip))
	for ii=1:Nsites
		phi_dist[ii,ip] = dum[ii]
	end
end
return phi_dist
end
################################################################################
function binder(psi,Nbasis,Xmat,Ymat)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis

include("operators.jl")


#Build O^2 operators#
Nsites = length(psi)
sites = siteinds(psi)

ampoX = AutoMPO()
ampoY = AutoMPO()

for i=1:Nsites
	ampoX += 1.0,"X",i
	ampoY += 1.0,"Y",i
end

Mx = MPO(ampoX,sites)
My = MPO(ampoY,sites)

#M2*Psi#
MxPsi = replaceprime(contract(Mx, psi), 2 => 1)#apply(Mx,psi)
Mx2Psi = replaceprime(contract(Mx, MxPsi), 2 => 1)#apply(Mx,MxPsi)
MyPsi = replaceprime(contract(My, psi), 2 => 1)#apply(My,psi)
My2Psi = replaceprime(contract(My, MyPsi), 2 => 1)#apply(My,MyPsi)

#Expectation values#
Mx2 = real(inner(psi,Mx2Psi))
Mx4 = real(inner(Mx2Psi,Mx2Psi))

My2 = real(inner(psi,My2Psi))
My4 = real(inner(My2Psi,My2Psi))

binderX = 1.0-Mx4/(3.0*Mx2^2)
binderY = 1.0-My4/(3.0*My2^2)

return Mx2,My2,binderX,binderY
end
###############################################################################################
function binder2(psi,Nbasis,Xmat,Xmat2,Xmat3,Xmat4,Ymat,Ymat2,Ymat3,Ymat4)

############################################################################################
##Calculation of binder parameter by using the multinomial theorem to decompose the powers##
############################################################################################

global X = Xmat
global X2 = Xmat2
global X3 = Xmat3
global X4 = Xmat4
global Y = Ymat
global Y2 = Ymat2
global Y3 = Ymat3
global Y4 = Ymat4
global Nspec = Nbasis

include("operators.jl")

Nsites = length(psi)

#Define list with operator pairings#
list1,list2,list3,list4 = index_list(Nsites)

######################
##<M^4>=<(sum x)^4>##
######################
Mx4 = 0.0
My4 = 0.0

#sum x^4#
x4 = expect(psi,"X4")
Mx4 += sum(x4)
y4 = expect(psi,"Y4")
My4 += sum(y4)

#sum xa^3*xb#
x3x = correlation_matrix(psi,"X3","X")
y3y = correlation_matrix(psi,"Y3","Y")
for ind in list1
	Mx4 += 4*x3x[ind[1],ind[2]]
	My4 += 4*y3y[ind[1],ind[2]]
end

#sum xa^2*xb^2#
x2x2 = correlation_matrix(psi,"X2","X2")
y2y2 = correlation_matrix(psi,"Y2","Y2")
for ind in list2
	Mx4 += 6*x2x2[ind[1],ind[2]]
	My4 += 6*y2y2[ind[1],ind[2]]
end

#sum xa^2*xb*xc#
x2xx,y2yy = threesite_correlation_x2xx(psi,Nbasis,Xmat,Xmat2,Ymat,Ymat2,list3)
Mx4 += 12*sum(x2xx)
My4 += 12*sum(y2yy)

#sum xa*xb*xc*xd#
xxxx,yyyy = foursite_correlation_xxxx(psi,Nbasis,Xmat,Ymat,list4)
Mx4 += 24*sum(xxxx)
My4 += 24*sum(yyyy)

#####################
##<M^2>=<(sum x)^2>##
#####################
Mx2 = 0.0
My2 = 0.0

#sum xa^2#
x2 = expect(psi,"X2")
Mx2 += sum(x2)
y2 = expect(psi,"Y2")
My2 += sum(y2)

#sum xa*xb#
xx = correlation_matrix(psi,"X","X")
yy = correlation_matrix(psi,"Y","Y")
for ind in list2
	Mx2 += 2*xx[ind[1],ind[2]]
	My2 += 2*yy[ind[1],ind[2]]
end

binderX=1.0-Mx4/(3.0*Mx2^2)
binderY=1.0-My4/(3.0*My2^2)

return binderX,binderY
end
################################################################################
function index_list(N)

#Define lists with term indices#
list1 = Tuple{Int,Int}[]
list2 = Tuple{Int,Int}[]
list3 = Tuple{Int,Int,Int}[]
list4 = Tuple{Int,Int,Int,Int}[]

for ii=1:N
	for jj=1:N
		if jj != ii
			push!(list1,(ii,jj))
		end
		if ii < jj
			push!(list2,(ii,jj))
		end
		for k=1:N
			if ii != jj && ii != k && jj != k && jj < k
				push!(list3,(ii,jj,k))
			end	
			for m=1:N
				if ii < jj < k < m
					push!(list4,(ii,jj,k,m))
				end	
			end
		end
	end
end

return list1,list2,list3,list4
end
################################################################################
function threesite_correlation_x2xx(psi,Nbasis,Xmat,Xmat2,Ymat,Ymat2,list3)

global X = Xmat
global X2 = Xmat2
global Y = Ymat
global Y2 = Ymat2
global Nspec = Nbasis

include("operators.jl")
Nsites = length(psi)

#sum xa^2*xb*xc#

##Add position of square term and sort list to reduce number of orthogonalizations##
list_x2xx = []
m=1
for ind in list3
	k=ind[1]
	l=sort(list3[m])
	push!(list_x2xx,(l[1],l[2],l[3],k))
	m+=1
end
sort!(list_x2xx,by = x -> x[1]) 


k1=0
k2=0
k3=0
t1=""
t2=""
t3=""
type = Array{String}(undef,3)
sites = siteinds(psi)
corr3_x = []
corr3_y = []
for ind in list_x2xx
	i1 = ind[1]
	i2 = ind[2]
	i3 = ind[3]
	pos = ind[4]
	type .="x"
	if i1 == pos
		type[1] = "x2"
	elseif i2 == pos
		type[2] = "x2"
	elseif i3 == pos
		type[3] = "x2"
	end

	###################################
	##Act with operator on first site##
	###################################
	#Define orthogonality center (only if first site changes)#
	if i1 != k1 || type[1] != t1 ##change in first site or operator type
		if i1 != k1 
			orthogonalize!(psi, i1) 
		end		
		#Determine which operator sits on site i1#
		if type[1] == "x2" 
			O1x = op("X2",sites[i1])
			O1y = op("Y2",sites[i1])
		else
			O1x = op("X",sites[i1])
			O1y = op("Y",sites[i1])
		end
		if i1 > 1
			rind = commonind(psi[i1],psi[i1+1])
			global op1x = psi[i1]*O1x*dag(prime(prime(psi[i1],"Site"),rind))
			global op1y = psi[i1]*O1y*dag(prime(prime(psi[i1],"Site"),rind))
		else
			global op1x = psi[i1]*O1x*dag(prime(psi[i1]))
			global op1y = psi[i1]*O1y*dag(prime(psi[i1]))
		end
	end

	########################################
	##Link region between operator 1 and 2##
	########################################
	if i1 != k1 || i2 != k2 || t1 != type[1]#first, second site or first operator type changed#
		global lr12x = copy(op1x)
		global lr12y = copy(op1y)
		for ii=(i1+1):(i2-1)
			global lr12x *= psi[ii]*dag(prime(psi[ii],"Link"))
			global lr12y *= psi[ii]*dag(prime(psi[ii],"Link"))
		end
	end
	####################################
	##Act with operator on second site##
	####################################
	if i1 != k1 || i2 != k2 || type[2] != t2 #new orthogonalization or second site/type changed#
		#Determine which operator sits on site i2#
		if type[2] == "x2" 
			O2x = op("X2",sites[i2])
			O2y = op("Y2",sites[i2])
		else
			O2x = op("X",sites[i2])
			O2y = op("Y",sites[i2])
		end
		tmp= psi[i2]*O2x*dag(prime(psi[i2]))
		global op2x = lr12x*tmp
		tmp= psi[i2]*O2y*dag(prime(psi[i2]))
		global op2y = lr12y*tmp
	end

	########################################
	##Link region between operator 2 and 3##
	########################################
	if i1 != k1 || i2 != k2 || i3 != k3 || t2 != type[2] #orthogonalization, second or third site changed#
		global lr23x = copy(op2x)
		global lr23y = copy(op2y)
		for ii=(i2+1):(i3-1)
			global lr23x *= psi[ii]*dag(prime(psi[ii],"Link"))
			global lr23y *= psi[ii]*dag(prime(psi[ii],"Link"))
		end
		k2=i2
	end

	###################################
	##Act with operator on third site##
	###################################
	if i1 != k1 || i3 != k3 || type[3] != t3 #orthogonalization, third site/type changed
		k3 = i3
		#Determine which operator sits on site i3#
		if type[3] == "x2" 
			O3x = op("X2",sites[i3])
			O3y = op("Y2",sites[i3])
		else
			O3x = op("X",sites[i3])
			O3y = op("Y",sites[i3])
		end
		if i3 < Nsites
			lind = commonind(psi[i3-1],psi[i3])
			global op3x = psi[i3]*O3x*dag(prime(prime(psi[i3],"Site"),lind))
			global op3y = psi[i3]*O3y*dag(prime(prime(psi[i3],"Site"),lind))
		else
			global op3x = psi[i3]*O3x*dag(prime(psi[i3]))
			global op3y = psi[i3]*O3y*dag(prime(psi[i3]))
		end
	end
	push!(corr3_x,scalar(lr23x*op3x))
	push!(corr3_y,scalar(lr23y*op3y))
	t1 = type[1]
	t2 = type[2]
	t3 = type[3]
	k1=i1
end
return corr3_x,corr3_y
end
################################################################################
function foursite_correlation_xxxx(psi,Nbasis,Xmat,Ymat,list4)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis
	
include("operators.jl")
Nsites = length(psi)

#sum xa*xb*xc*xd#
k1=0
k2=0
k3=0
k4=0
corr4_x=[]
corr4_y=[]
sites = siteinds(psi)
for ind in list4
	i1 = ind[1]
	i2 = ind[2]
	i3 = ind[3]
	i4 = ind[4]

	###################################
	##Act with operator on first site##
	###################################
	#Define orthogonality center (only if first site changes)#
	if i1 != k1 ##change in first site 
		orthogonalize!(psi, i1) 	
		O1x = op("X",sites[i1])
		O1y = op("Y",sites[i1])
		if i1 > 1
			rind = commonind(psi[i1],psi[i1+1])
			global op1x = psi[i1]*O1x*dag(prime(prime(psi[i1],"Site"),rind))
			global op1y = psi[i1]*O1y*dag(prime(prime(psi[i1],"Site"),rind))
		else
			global op1x = psi[i1]*O1x*dag(prime(psi[i1]))
			global op1y = psi[i1]*O1y*dag(prime(psi[i1]))
		end
	end

	########################################
	##Link region between operator 1 and 2##
	########################################
	if i1 != k1 || i2 != k2 #first or second site changed#
		global lr12x = copy(op1x)
		global lr12y = copy(op1y)
			for ii=(i1+1):(i2-1)
				global lr12x *= psi[ii]*dag(prime(psi[ii],"Link"))
				global lr12y *= psi[ii]*dag(prime(psi[ii],"Link"))
			end
	end

	####################################
	##Act with operator on second site##
	####################################
	if i1 != k1 || i2 != k2 #new orthogonalization or second site changed#
		O2x = op("X",sites[i2])
		O2y = op("Y",sites[i2])
		tmp= psi[i2]*O2x*dag(prime(psi[i2]))
		global op2x = lr12x*tmp
		tmp= psi[i2]*O2y*dag(prime(psi[i2]))
		global op2y = lr12y*tmp
	end

	########################################
	##Link region between operator 2 and 3##
	########################################
	if i1 != k1 || i2 != k2 || i3 != k3 #orthogonalization, second or third site changed#
		global lr23x = copy(op2x)
		global lr23y = copy(op2y)
		for ii=(i2+1):(i3-1)
			global lr23x *= psi[ii]*dag(prime(psi[ii],"Link"))
			global lr23y *= psi[ii]*dag(prime(psi[ii],"Link"))
		end
	end	

	####################################
	##Act with operator on third site##
	####################################
	if i1 != k1 || i2 != k2 || i3 != k3 #new orthogonalization, second or third site changed#
		O3x = op("X",sites[i3])
		O3y = op("Y",sites[i3])
		tmp= psi[i3]*O3x*dag(prime(psi[i3]))
		global op3x = lr23x*tmp
		tmp= psi[i3]*O3y*dag(prime(psi[i3]))
		global op3y = lr23y*tmp
	end

	########################################
	##Link region between operator 3 and 4##
	########################################
	if i1 != k1 || i2 != k2 || i3 != k3 || i4 != k4 #orthogonalization, third or fourth site changed#
		global lr34x = copy(op3x)
		global lr34y = copy(op3y)
		for ii=(i3+1):(i4-1)
			global lr34x *= psi[ii]*dag(prime(psi[ii],"Link"))
			global lr34y *= psi[ii]*dag(prime(psi[ii],"Link"))
		end
		k3=i3
	end	

	####################################
	##Act with operator on fourth site##
	####################################
	if i1 != k1 || i4 != k4 #orthogonalization, third site changed
		k4 = i4
		O4x = op("X",sites[i4])
		O4y = op("Y",sites[i4])
		if i4 < Nsites
			lind = commonind(psi[i4-1],psi[i4])
			global op4x = psi[i4]*O4x*dag(prime(prime(psi[i4],"Site"),lind))
			global op4y = psi[i4]*O4y*dag(prime(prime(psi[i4],"Site"),lind))
		else
			global op4x = psi[i4]*O4x*dag(prime(psi[i4]))
			global op4y = psi[i4]*O4y*dag(prime(psi[i4]))
		end
	end
	push!(corr4_x,scalar(lr34x*op4x))
	push!(corr4_y,scalar(lr34y*op4y))
	k1=i1	
	k2=i2
end


return corr4_x,corr4_y
end
################################################################################
end
