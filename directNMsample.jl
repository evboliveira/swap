using ITensors,ITensorMPS
using Printf
using LinearAlgebra
using Statistics
import StatsBase as SB
using HypothesisTests
import Distributions as Dist
using FFTW
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
#### Int Strength g values ####
listg = append!(
    [g for g = 0.0 : 0.1 : 0.4],
    [g for g = 0.41 : 0.01 : 0.60],
    [g for g = 0.7 : 0.1 : 5.0])
Ng = length(listg)
################################
mmax = 10
mbond = Nsites ÷ 2 #for the bipartite entanglement
pairs = "nearest" # nearest,allpairs
evod = "all"	# all,dvr,all_real

#### RESULTS PATH ###

# ### Locally
# res_path = "./results/N$Nsites/"

### Compute Canada
res_path = "/home/evbdeoli/scratch/results/swap/28may2026/N$Nsites/"
######################

#Define output files#
create_file(res_path*"directNMswap.txt")

#### SAMPLING PROCEDURE #####
Nsamples = 10000

for ig = 0:length(listg)-1
	let
		g= listg[ig+1]
		println(string("#### Sampling g=",string(round(g,digits=3))))
		mps_out=h5open(string(res_path*"psi0/psi0_g",string(round(g,digits=3))),"r")
		psi=read(mps_out,"MPS",MPS)
		close(mps_out)
        
		sites = siteinds(psi)

		directNMswap = zeros(Float64,Nsamples)

		for conf=1:Nsamples			
			rep1=sample(psi) ## replica 1
			rep2=sample(psi) ## replica 2

			#### direct NM entanglement entropy

			exponent=0
			for i=1:Nsites
				for k=1:Nsites
					if k != i
						exponent += (rep1[i]-(mmax+1))*(rep1[k]-(mmax+1)) + (rep2[i]-(mmax+1))*(rep2[k]-(mmax+1)) - 2*(rep1[i]-(mmax+1))*(rep2[k]-(mmax+1))
					end
				end
			end
			directNMswap[conf] = exp(exponent)	
		end

		directNMswap_avg = mean(directNMswap)
		directNM_S2 = -log(directNMswap_avg)
		directNMerror = sqrt(var(directNMswap)/Nsamples)/directNMswap_avg
        
		write_output_error(res_path*"directNMswap.txt",g,directNM_S2,directNMerror)
	end
end