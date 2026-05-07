using ITensors,ITensorMPS
using Printf
using LinearAlgebra
using Statistics
import StatsBase as SB
using HypothesisTests
import Distributions as Dist
using FFTW
using Plots
using HDF5
push!(LOAD_PATH,pwd())
push!(LOAD_PATH,"../")
using matrices
using states
using expectations
using distributions
using create_operators 
using dvr
using sign

function create_file(path)
	f=open(path,"w")
	if evod == "dvr"
		println(f,"#Ngrid= ",mmax)
		println(f,"#DVR basis")
	else
		println(f,"#mmax= ",mmax)
		println(f,"#m-states: ",evod," states (FBR)")
	end
	println(f,"#Nr. of sites: ",Nsites)
	println(f)
	close(f)
end

function write_output(path,g,observable)
        text=' '
        for b=1:length(observable)
                text*=string(observable[b],' ')
        end
	f=open(path,"a")
	println(f,round(g,digits=4)," ",text)
	close(f)
end

function write_output_error(path,g,observable,error)
        text=' '
        for b=1:length(observable)
                text*=string(observable[b],' ')
        end
	f=open(path,"a")
	println(f,round(g,digits=4)," ",text," ",error)
	close(f)
end

function write_text(path,g,txt)
	f=open(path,"a")
	println(f,round(g,digits=4)," ",txt)
	close(f)
end

function transform_realbasis(A,U)
	tmp = BLAS.gemm('C','N', U,A)
	B = BLAS.gemm('N','N', tmp,U)

	return B
end

function trans_matrix(mmax)
	Utrans = zeros(ComplexF64,(2*mmax+1,2*mmax+1))

	for m=1:mmax
		Utrans[m,m] = 1.0/sqrt(2.0)
		Utrans[m,2*mmax+2-m] = 1.0im/sqrt(2.0)
	end
	Utrans[mmax+1,mmax+1] = 1.0
	for m=1:mmax
		Utrans[2*mmax+2-m,m] = 1.0/sqrt(2.0)
		Utrans[2*mmax+2-m,2*mmax+2-m] = -1.0im/sqrt(2.0)
	end
	
	return Utrans
end


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

f=open("log","w")

println(f,"#######################################")
println(f,"###########Basis information###########")
println(f,"#######################################")
println(f,"mmax= ",mmax)
println(f,"Number of sites: ",Nsites)
println(f,"Dimension of local Hilbert space: ",2*mmax+1)
println(f,"V6 Strength ",V6strength)

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
	println(f,"all m-states are considered")
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
	println(f,"all m-states are considered")
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
	println(f,"DVR-basis is used")
end

#Calculate higher powers of X#
global X2 = BLAS.gemm('N','N', X,X)
global X3 = BLAS.gemm('N','N', X2,X)
global X4 = BLAS.gemm('N','N', X3,X)
global Y2 = BLAS.gemm('N','N', Y,Y)
global Y3 = BLAS.gemm('N','N', Y2,Y)
global Y4 = BLAS.gemm('N','N', Y3,Y)

if pairs == "nearest"
	println(f,"only nearest-neighbour interactions")
elseif pairs == "allpairs"
	println(f,"all interactions")
end
#Determine number of interaction pairs per starting site#
Nsecond = zeros(Int64,(Nsites-1))
for i=1:Nsites-1
        if pairs == "nearest"
                Nsecond[i]=i+1
        elseif pairs == "allpairs"
                Nsecond[i]=Nsites
        end
end


println(f)
println(f,"#################################################################################")
println(f,"##################################")
println(f,"####Calculate free rotor chain####")
println(f,"##################################")
println(f)
close(f)

#Define output files#
create_file("L.txt")
create_file("entropy_swap.txt")
create_file("entropy_swap0.txt")
create_file("swap0.txt")



listg=[]

for ii=0:4
	push!(listg,0.1*ii)
end

for ii=1:20
	push!(listg,0.4+0.01*ii)
end

for ii=1:24
	push!(listg,0.6+0.1*ii)
end

#listg=[0.5,0.7,1.1,1.5,2.0]
#for ig = 0:Ng-1
f=open("sample","w")

fast=0
for ig = 0:length(listg)-1
	let
		g= listg[ig+1]
		mps_out=h5open(string("psi0_g",string(round(g,digits=3))),"r")
		psi=read(mps_out,"MPS",MPS)
		close(mps_out)

		data_out=open(string("L",string(round(g,digits=3))),"w")
		sites = siteinds(psi)

		Na=mbond
		#psi = orthogonalize(psi, 1)
		Nsamples=4000


		sample1_a=zeros(Int64,(Na))
		sample1_b=zeros(Int64,(Na))
		sample2_a=zeros(Int64,(Na))
		sample2_b=zeros(Int64,(Na))

		qsample1_a=zeros(Float64,(1))
		qsample1_b=zeros(Float64,(Nsites-1))
		qsample2_a=zeros(Float64,(1))
		qsample2_b=zeros(Float64,(Nsites-1))

		sample_swap1=zeros(Int64,(Nsites))
		sample_swap2=zeros(Int64,(Nsites))

		qsample_swap1=zeros(Float64,(Nsites))
		qsample_swap2=zeros(Float64,(Nsites))

		inv_qsample_swap1=zeros(Int64,(Nsites))
		inv_qsample_swap2=zeros(Int64,(Nsites))

		swap=zeros(Float64,Nsamples)
		swap0=zeros(Float64,Nsamples)
		totalL=zeros(Float64,Nsamples)
		totalLbin=zeros(Int64,Nsamples)

		nbin=(2*mmax+1)*Nsites

		for conf=1:Nsamples
			s1=sample(psi)
			L=0
			Lbin=0
			for i=1:Nsites
				L+=s1[i]-mmax-1
				Lbin+=s1[i]
			end
			totalL[conf]=L
			totalLbin[conf]=Lbin/2
			s2=sample(psi)

			for i=1:Na
				sample1_a[i]= s1[i]
				sample2_a[i]= s2[i]
				sample1_b[i]= s1[i+Na]
				sample2_b[i]= s2[i+Na]
			end

			qs1=dct(s1)
			qs2=dct(s2)
			qsample1_a[1]= qs1[1]
			qsample2_a[1]= qs2[1]
			for i=2:Nsites
				qsample1_b[i-1]= qs1[i]
				qsample2_b[i-1]= qs2[i]
			end

			for i=1:Na
				sample_swap1[i]=sample2_a[i]
				sample_swap1[i+Na]=sample1_b[i]
				sample_swap2[i]=sample1_a[i]
				sample_swap2[i+Na]=sample2_b[i]
			end
			qsample_swap1[1]=qsample2_a[1]
			qsample_swap2[1]=qsample1_a[1]
			for i=2:Nsites
				qsample_swap1[i]=qsample1_b[i-1]
				qsample_swap2[i]=qsample2_b[i-1]
			end

			qs1=qsample_swap1[:]
			qs2=qsample_swap2[:]
			st1=idct(qs1)
			st2=idct(qs2)
			for i=1:Nsites
				if st1[i]>(2*mmax+1)
					println("basis state out of range",st1[i])
					st1[i]=2*mmax+1
				end
				if st1[i]<1
					println("basis state out of range",st1[i])
					st1[i]=1
				end
				if st2[i]>(2*mmax+1)
					println("basis state out of range",st2[i])
					st2[i]=2*mmax+1
				end
				if st2[i]<1
					println("basis state out of range",st2[i])
					st2[i]=1
				end
				inv_qsample_swap1[i]=round(Int,st1[i])
				inv_qsample_swap2[i]=round(Int,st2[i])
				#println(inv_qsample_swap1[conf,i]-mmax)
				#println(inv_qsample_swap2[conf,i]-mmax)
			end

			sw1=sample_swap1[:]
			if fast==0
				V = ITensor(1.)
				for j=1:Nsites
					V *= (psi[j]*state(sites[j],sw1[j]))
				end
				wf1 = scalar(V)
			else
				wf1=inner(MPS(sites,sw1),psi)
			end

			sw2=sample_swap2[:]
			if fast==0
				V = ITensor(1.)
				for j=1:Nsites
					V *= (psi[j]*state(sites[j],sw2[j]))
				end
				wf2 = scalar(V)
			else
				wf2=inner(MPS(sites,sw2),psi)
			end

			if fast==0
				V = ITensor(1.)
				for j=1:Nsites
					V *= (psi[j]*state(sites[j],s1[j]))
				end
				wf3 = scalar(V)
			else
				wf3=inner(MPS(sites,s1),psi)
			end
			
			if fast==0
				V = ITensor(1.)
				for j=1:Nsites
					V *= (psi[j]*state(sites[j],s2[j]))
				end
				wf4 = scalar(V)
			else
				wf4=inner(MPS(sites,s2),psi)
			end

			swap[conf]=real(wf1*wf2/wf3/wf4)
			#println(f,swap[conf]," ",wf1," ",wf2," ",wf3," ",wf4);flush(f)

			# total m
			sw1=inv_qsample_swap1[:]
			V = ITensor(1.)
			for j=1:Nsites
				V *= (psi[j]*state(sites[j],sw1[j]))
			end
			wf1 = scalar(V)
			sw2=inv_qsample_swap2[:]
			V = ITensor(1.)
			for j=1:Nsites
				V *= (psi[j]*state(sites[j],sw2[j]))
			end
			wf2 = scalar(V)

			value=real(wf1*wf2/wf3/wf4)
			if value >0 && value <= 10
				swap0[conf]=value
			end
			#println(f,swap0[conf]," ",value," ",real(wf1)," ",real(wf2)," ",real(wf3)," ",real(wf4));flush(f)
			#println(f,swap0[conf]," ",s1," ",s2," ",s3," ",s4);flush(f)

		end
		swap_avg = mean(swap)
		swap0_avg = mean(swap0)
		S2 = -log(swap_avg)
		S20 = -log(swap0_avg)
		error = sqrt(var(swap)/Nsamples)/swap_avg
		error0 = sqrt(var(swap0)/Nsamples)/swap0_avg

		# Define theoretical Gaussian distribution (mean=0, std=1)
		d = Dist.Normal(mean(totalLbin), sqrt(var(totalLbin)))

		for idx=1:Nsamples
			println(data_out,totalL[idx])
		end

		close(data_out)
		# fitted_dist = Dist.fit_mle(Dist.Binomial,nbin,totalLbin)
		# println(fitted_dist)

		# # Sample stats
		# println("Sample Mean: ", mean(totalLbin))
		# println("Sample Var:  ", var(totalLbin))

		# # Distribution stats
		# println("Theory Mean: ", mean(fitted_dist))
		# println("Theory Var:  ", var(fitted_dist))

		# Perform the K-S test
		ks_test = ExactOneSampleKSTest(totalLbin, d)

		#histogram(totalL)

		observable=[mean(totalL),sqrt(var(totalL)/Nsamples),var(totalL),SB.skewness(totalL),SB.kurtosis(totalL),pvalue(ks_test)]
		write_output("L.txt",g,observable)
		write_output_error("entropy_swap.txt",g,S2,error)
		write_output_error("entropy_swap0.txt",g,S20,error0)
		write_output_error("swap0.txt",g,swap0_avg,sqrt(var(swap0)/Nsamples))

	end
end
