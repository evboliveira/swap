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
push!(LOAD_PATH,"src/")
# push!(LOAD_PATH,pwd())
# push!(LOAD_PATH,"../")
using matrices
using states
using expectations
using distributions
using create_operators 
using dvr
using sign

#### AUXILIARY FUNCTIONS #### # !!! Move to the module in src/aux_suncs !!!
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
########################################################

#### INPUT PARAMETERS ####
# System
Nsites = parse(Int, ARGS[1]) ## as an external argument
listg = append!(
    [g for g = 0.0 : 0.1 : 0.4],
    [g for g = 0.41 : 0.01 : 0.60],
    [g for g = 0.7 : 0.1 : 3.0])
Ng = length(listg)
mmax = 5
mbond = Nsites ÷ 2 #for the bipartite entanglement
pairs = "nearest" # nearest,allpairs
evod = "all"	# all,dvr,all_real

#### RESULTS PATH ###

### Locally
res_path = "./results/N$Nsites/"

### Compute Canada
# res_path = "/home/evbdeoli/scratch/results/swap/N$Nsites/"
######################

#Define output files#
create_file(res_path*"L.txt")
create_file(res_path*"entropy_swap.txt")
create_file(res_path*"NMentropy_swap.txt")
create_file(res_path*"NMswap.txt")

#### Int Strength g values ####
listg=[]
# listg = [g for g = 0.1 : 0.1 : 2.0]
listg = append!(
    [g for g = 0.0 : 0.1 : 0.4],
    [g for g = 0.41 : 0.01 : 0.60],
    [g for g = 0.7 : 0.1 : 3.0])
Ng = length(listg)
################################

#### SAMPLING PROCEDURE #####
Nsamples = 4000
Na=mbond

fast = true #Fast sampling?

for ig = 0:length(listg)-1
	let
		g= listg[ig+1]
		println(string("#### Sampling g=",string(round(g,digits=3))))
		mps_out=h5open(string(res_path*"psi0/psi0_g",string(round(g,digits=3))),"r")
		psi=read(mps_out,"MPS",MPS)
		close(mps_out)

		data_out=open(string(res_path*"L",string(round(g,digits=3))),"w")
		sites = siteinds(psi)
		#psi = orthogonalize(psi, 1)

		## Lattice Bipartite samples ##
		latt_swap1=zeros(Int64,(Nsites))
		latt_swap2=zeros(Int64,(Nsites))

		swap = zeros(Float64,Nsamples)
		####

		## Lattice NM samples (0th mode vs the rest Rth)##
		NMlatt_swap1=zeros(Float64,(Nsites))
		NMlatt_swap2=zeros(Float64,(Nsites))

		inv_NMlatt_swap1=zeros(Int64,(Nsites))
		inv_NMlatt_swap2=zeros(Int64,(Nsites))
		
		NMswap = zeros(Float64,Nsamples)
		####

		## Total Ang. Mom. stats. ##
		totalL=zeros(Float64,Nsamples)
		totalLbin=zeros(Int64,Nsamples)		
		nbin=(2*mmax+1)*Nsites
		####

		for conf=1:Nsamples			
			rep1=sample(psi) ## replica 1
			rep2=sample(psi) ## replica 2

			## Total Ang Mom Stats ##
			L = sum(rep1 .-(mmax+1))
			Lbin = sum(rep1)
			totalL[conf]=L
			totalLbin[conf]=Lbin ÷ 2
			#####

			## Lattice Bipartite sampling ##
			latt_swap1[Na+1:Nsites] = rep1[Na+1:Nsites] ## 1B --> 1B
			latt_swap1[1:Na] = rep2[1:Na] 			   ## 1A --> 2A
			latt_swap2[1:Na] = rep1[1:Na]              ## 2A --> 1A
			latt_swap2[Na+1:Nsites] = rep2[Na+1:Nsites] ## 2B --> 2B

			if fast
				V1 = ITensor(1.)
				V2 = ITensor(1.)
				V3 = ITensor(1.)
				V4 = ITensor(1.)
				for j=1:Nsites
					V1 *= (psi[j]*state(sites[j],latt_swap1[j]))
					V2 *= (psi[j]*state(sites[j],latt_swap2[j]))
					V3 *= (psi[j]*state(sites[j],rep1[j]))
					V4 *= (psi[j]*state(sites[j],rep2[j]))
				end
				wf1 = scalar(V1)
				wf2 = scalar(V2)
				wf3 = scalar(V3)
				wf4 = scalar(V4)
			else
				wf1=inner(MPS(sites,latt_swap1),psi)
				wf2=inner(MPS(sites,latt_swap2),psi)
				wf3=inner(MPS(sites,rep1),psi)
				wf4=inner(MPS(sites,rep2),psi)
			end

			swap[conf]=real(wf1*wf2/wf3/wf4)
			################################


			## NM lattice sampling ##

			## discrete cosine transforms ##
			NMrep1 = dct(rep1) ## NM replica 1
			NMrep2 = dct(rep2) ## NM replica 2
			##

			NMlatt_swap1[2:Nsites] = NMrep1[2:Nsites]  ## 1_Rth --> 1_Rth
			NMlatt_swap1[1] = NMrep2[1] 			   ## 1_0th --> 2_0th
			NMlatt_swap2[1] = NMrep1[1]                ## 2_0th --> 1_0th
			NMlatt_swap2[2:Nsites] = NMrep2[2:Nsites]  ## 2_Rth --> 2_Rth

			## Inverse discrete cosine transforms ##
			NMaux1 = idct(NMlatt_swap1)
			NMaux2 = idct(NMlatt_swap2)

			## Regularizing the basis states:
			for i=1:Nsites
				if NMaux1[i] > (2*mmax+1)
					println("basis state out of range",NMaux1[i])
					NMaux1[i]=2*mmax+1
				elseif NMaux1[i] < 1
					println("basis state out of range",NMaux1[i])
					NMaux1[i]=1
				end
				if NMaux2[i] > (2*mmax+1)
					println("basis state out of range",NMaux2[i])
					NMaux2[i]=2*mmax+1
				elseif NMaux2[i] < 1
					println("basis state out of range",NMaux2[i])
					NMaux2[i]=1
				end
				inv_NMlatt_swap1[i]=round(Int,NMaux1[i])
				inv_NMlatt_swap2[i]=round(Int,NMaux2[i])
			end

			V1 = ITensor(1.)
			V2 = ITensor(1.)
			for j=1:Nsites
				V1 *= (psi[j]*state(sites[j],inv_NMlatt_swap1[j]))
				V2 *= (psi[j]*state(sites[j],inv_NMlatt_swap2[j]))
			end
			wf1 = scalar(V1)
			wf2 = scalar(V2)

			value = real(wf1*wf2/wf3/wf4)
			if value >0 && value <= 10
				NMswap[conf]=value
			end
			
		end
		swap_avg = mean(swap)
		NMswap_avg = mean(NMswap)
		S2 = -log(swap_avg)
		NM_S2 = -log(NMswap_avg)
		error = sqrt(var(swap)/Nsamples)/swap_avg
		NMerror = sqrt(var(NMswap)/Nsamples)/NMswap_avg

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

		observable=[mean(totalL),sqrt(var(totalL)/Nsamples),var(totalL),SB.skewness(totalL),SB.kurtosis(totalL),pvalue(ks_test)]
		write_output(res_path*"L.txt",g,observable)
		write_output_error(res_path*"entropy_swap.txt",g,S2,error)
		write_output_error(res_path*"NMentropy_swap.txt",g,NM_S2,NMerror)
		write_output_error(res_path*"NMswap.txt",g,NMswap_avg,sqrt(var(NMswap)/Nsamples))

	end
end
