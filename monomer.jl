using Printf
using LinearAlgebra
push!(LOAD_PATH,pwd())
using matrices
using states
using expectations
using distributions
using create_operators 
using dvr
using sign




f=open("input.txt")
lines=readlines(f)
close(f)
lsplit=split(lines[2])
mmax=parse(Int64, lsplit[2])
lsplit=split(lines[5])
Nsites=parse(Int64, lsplit[2])
lsplit=split(lines[6])
Nbonds=parse(Int64, lsplit[2])
lsplit=split(lines[7])
Nsweep=parse(Int64, lsplit[2])
lsplit=split(lines[8])
e_cutoff=parse(Float64, lsplit[2])
lsplit=split(lines[9])
SVD_error=parse(Float64, lsplit[2])
lsplit=split(lines[12])
gstart=parse(Float64, lsplit[2])
lsplit=split(lines[13])
delta_g=parse(Float64, lsplit[2])
lsplit=split(lines[14])
Ng=parse(Int64, lsplit[2])
lsplit=split(lines[15])
mbond=parse(Int64, lsplit[2])
lsplit=split(lines[18])
pairs=lsplit[2]
lsplit=split(lines[21])
evod=lsplit[2]
lsplit=split(lines[24])
angle=parse(Float64, lsplit[2])
angle=angle*pi/180.0
lsplit=split(lines[25])
Estrength=parse(Float64, lsplit[2])
lsplit=split(lines[28])
Nstates=parse(Int64, lsplit[2])
lsplit=split(lines[31])
V6strength=parse(Float64, lsplit[2])


listV6=[]

for ii=0:10
	push!(listV6,0.1*ii)
end

for ii=1:99
	push!(listV6,1+ii)
end

for ii=1:5
	push!(listV6,100+ii*100)
end

for V6strength in listV6
	fV6=open("V6_mono"*string(V6strength)*".dat","w")

	fV6ev=open("V6_mono"*string(V6strength)*"ev.dat","w")

	Nm=30
	Elist=zeros((2*Nm+1,2*Nm+1))
	for mmax=1:Nm

		d=2*mmax+1

		#Calculate kinetic matrix and x operator#
		Ttmp = kinetic(mmax)
		Xtmp = Xoperator(mmax)
		Ytmp = Yoperator(mmax)
		Lztmp = lz_operator(mmax)
		V6tmp = V6_operator(mmax)

		#Define basis#
		if evod == "all"
			Nphi=90
			global phi = [ii*2.0*pi/Nphi for ii=1:Nphi] 
			Dtmp = distro_complex(mmax,phi) 
			global T = Ttmp
			global X = Xtmp
			global Y = Ytmp
			global D = Dtmp
			global Lz = Lztmp 
			global V6 = V6tmp 
			Nspec=size(T,1)
		elseif evod == "all_real"	
			
			#Real basis#
			Utrans = trans_matrix(mmax)
			Ttmp2 = real(transform_realbasis(Ttmp .+ 0.0im,Utrans))
			Xtmp2 = real(transform_realbasis(Xtmp .+ 0.0im,Utrans))
			Ytmp2 = real(transform_realbasis(1.0im*Ytmp,Utrans))
			Lztmp2 = imag(transform_realbasis(Lztmp .+ 0.0im,Utrans))
			
			Nphi=90
			global phi = [ii*2.0*pi/Nphi for ii=1:Nphi] 
			Dcomp = distro_complex(mmax,phi) 
			Dtmp = zeros(Float64,(2*mmax+1,2*mmax+1,Nphi))
			for ip=1:Nphi
				Dtmp[:,:,ip] = real(transform_realbasis(Dcomp[:,:,ip],Utrans))	
			end

			global T = Ttmp2
			global X = Xtmp2
			global Y = Ytmp2
			global Lz = Lztmp2
			global D = Dtmp
			Nspec=size(T,1)
		elseif evod == "dvr"
			tmp1,tmp2,tmp3,tmp4,tmp5 = exp_dvr(mmax)
			global T = tmp1
			global X = tmp2
			global Y = tmp3
			global Lz = imag(tmp4)
			global V6= tmp5
			Nphi=2*mmax+1
			global phi = [ii*2.0*pi/(2*mmax+1) for ii=1:(2*mmax+1)]
			global D = distro_dvr(mmax,phi) 
			Nspec=size(T,1)
		end

		#monomer diagonalization
		Hmono=T-0.5*V6strength*V6
		evals,evecs=eigen(Hmono)
		for i=1:d
			Elist[mmax,i]=evals[i]
		end

		if mmax==Nm
			npoints=100
			dphi=2*pi/npoints
			psi=zeros(ComplexF64,(npoints))
			for i=1:d
				for a =1:npoints
					phase=(-mmax+i-1)*(a-1)*dphi
					psi[a]+=evecs[i,1]*exp(1.0im*phase)
				end
			end
			for a =1:npoints
				println(fV6ev,(a-1)*dphi," ",abs(psi[a]))
			end
		end
	end
	for mmax=1:Nm
		d=2*mmax+1
		for i=1:d
			print(fV6,mmax," ",Elist[mmax,i]+0.5*V6strength)
		end
		println(fV6)
	end
	close(fV6)
	close(fV6ev)
end
