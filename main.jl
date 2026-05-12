using ITensors,ITensorMPS
using Printf
using LinearAlgebra
using Statistics
using FFTW
using HDF5
using Dates
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
using aux_funcs

#### AUXILIARY FUNCTIONS #### # !!! Move to the module in src/aux_suncs !!!
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
	println(f," ")
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
	log_println(f,round(g,digits=4)," ",text," ",error)
	close(f)
end

function write_text(path,g,txt)
	f=open(path,"a")
	log_println(f,round(g,digits=4)," ",txt)
	close(f)
end
########################################################

#### LOADING INPUT ####
f=open("input.txt")
lines=readlines(f)
close(f)
lsplit=split(lines[2])
mmax=parse(Int64, lsplit[2])
lsplit=split(lines[5])
# Nsites=parse(Int64, lsplit[2])
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
########

#### RESULTS PATH ####
res_path = "./results/"
time_stamp = Dates.format(now(),"dd-mm-yyyy_HH:MM:SS")
res_path *= time_stamp
mkdir(res_path)
res_path *= "/"
mkdir(res_path*"psi0")
######################


f=open(res_path*"log","w")
log_println(f,"#######################################")
log_println(f,"###########Basis information###########")
log_println(f,"#######################################")
log_println(f,"mmax= ",mmax)
log_println(f,"Number of sites: ",Nsites)
log_println(f,"Dimension of local Hilbert space: ",2*mmax+1)
log_println(f,"V6 Strength ",V6strength)

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
	log_println(f,"all m-states are considered")

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
	log_println(f,"all m-states are considered")
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
	log_println(f,"DVR-basis is used")
end

#Calculate higher powers of X#
global X2 = BLAS.gemm('N','N', X,X)
global X3 = BLAS.gemm('N','N', X2,X)
global X4 = BLAS.gemm('N','N', X3,X)
global Y2 = BLAS.gemm('N','N', Y,Y)
global Y3 = BLAS.gemm('N','N', Y2,Y)
global Y4 = BLAS.gemm('N','N', Y3,Y)

if pairs == "nearest"
	log_println(f,"only nearest-neighbour interactions")
elseif pairs == "allpairs"
	log_println(f,"all interactions")
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


println(f," ")
log_println(f,"#################################################################################")
log_println(f,"##################################")
log_println(f,"####Calculate free rotor chain####")
log_println(f,"##################################")
println(f," ")
close(f)

#Define output files#
create_file(res_path*"energy.txt")
create_file(res_path*"entropy_vN.txt")
create_file(res_path*"entropy_Renyi.txt")
# create_file(res_path*"entropy_swap.txt")
# create_file(res_path*"entropy_swap0.txt")
# create_file(res_path*"swap0.txt")
create_file(res_path*"mux.txt")
create_file(res_path*"muy.txt")
create_file(res_path*"xcorr.txt")
create_file(res_path*"ycorr.txt")
create_file(res_path*"corr.txt")
create_file(res_path*"mcorr.txt")
create_file(res_path*"binder_x.txt")
create_file(res_path*"binder_y.txt")
# create_file(res_path*"tau.txt")
# create_file(res_path*"tau_sample.txt")
# create_file(res_path*"cplane_sample.txt")
# create_file(res_path*"cplane.txt")
if Nstates > 0
	for i=1:Nstates+1
		create_file(res_path*"schmidt_values_"*string(i)*".txt")
		create_file(res_path*"phi_distribution_"*string(i)*".txt")
	end
else
	create_file(res_path*"schmidt_values.txt")
	create_file(res_path*"phi_distribution.txt")
end

operators(Nspec,Nphi,evod)

listg=[]
# listg = [g for g = 0.1 : 0.1 : 2.0]
listg = append!(
    [g for g = 0.0 : 0.1 : 0.4],
    [g for g = 0.41 : 0.01 : 0.60],
    [g for g = 0.7 : 0.1 : 3.0])
Ng = length(listg)

for ig = 0:Ng-1
let
	include("src/operators.jl")
	include("src/observer.jl")

	if evod == "dvr" || evod == "all_real"
		fac1 = 1.0
		fac2 = 1.0
	else 
		fac1 = -1.0
		fac2 = 1.0im
	end

	#Non-interacting rotors as initial guess#
	if ig == 0
		sites=siteinds("PlaRotor",Nsites)
		ampo0 = AutoMPO()
		#Define Hamiltonian as MPO#
		for j=1:Nsites
			ampo0 += 1.0,"T",j
			#Electric field#
			ampo0 += -cos(angle)*Estrength,"X",j
			ampo0 += -sin(angle)*Estrength*fac2,"Y",j
			ampo0 += (-0.5*V6strength),"V6",j
		end	
		H0 = MPO(ampo0,sites)
		#Define accuracy parameters#
		sweeps = Sweeps(Nsweep)
		#Set up initial state#
		global psi0 = randomMPS(sites,10)
		maxdim!(sweeps,10) # gradually increase states kept
		cutoff!(sweeps,SVD_error) # desired truncation error
		
		#Perform DMRG runs#
		obs = DemoObserver(e_cutoff)
		energy,psi = dmrg(H0,psi0,sweeps,observer=obs, outputlevel=0)
		
		global psi0 = psi
		maxbond=maxlinkdim(psi)
		f=open(res_path*"log","a")
		log_println(f,"Max. bond dimension: ",maxbond)
		println(f," ")
		@printf(f,"Final energy = %.8f \n",energy)
		println(f," ")
		log_println(f,"Initial state calculated")
		log_println(f,"###############################################################################")
		println(f," ")
		close(f)
	end	
	#g = gstart + ig*delta_g
	g= listg[ig+1]
	f=open(res_path*"log","a")
	log_println(f,"##################################")
	log_println(f,"########g= ",g," ########")
	log_println(f,"##################################")
	println(f," ")
	log_println(f,"####DMRG calculation####")
	log_println(f,"Construct MPO")
	

	sites = siteinds(psi0)
	ampo = AutoMPO() 
	for i=1:Nsites-1
		ampo += 1.0,"T",i
		ampo += (-0.5*V6strength),"V6",i

		for j=i+1:Nsecond[i]
			c=g/((abs(j-i))^3)
			#y_iy_j#
			ampo += 1.0*c*fac1,"Y",i,"Y",j
			#2*x_ix_j#
			ampo += -2.0*c,"X",i,"X",j
		end
		#Electric field#
		ampo += -cos(angle)*Estrength,"X",i
		ampo += -sin(angle)*Estrength*fac2,"Y",i
	end
	ampo += 1.0,"T",Nsites
	#Electric field#
	ampo += -cos(angle)*Estrength,"X",Nsites
	ampo += -sin(angle)*Estrength*fac2,"Y",Nsites
	ampo += (-0.5*V6strength),"V6",Nsites


	H = MPO(ampo,sites)
	#Define accuracy parameters#
	sweeps = Sweeps(Nsweep)
	#Set up initial state#
	cutoff!(sweeps,SVD_error) # desired truncation error
	if ig == 0
		maxdim!(sweeps,10,20,30,40,Nbonds) # gradually increase states kept
	else	
		maxdim!(sweeps,maxlinkdim(psi0),Nbonds)
	end

	#Perform DMRG runs#
	log_println(f,"Start DMRG run")
	close(f)
	obs = DemoObserver(e_cutoff)
	energy,psi = dmrg(H,psi0,sweeps,observer=obs, outputlevel=0)

		

	#################################################################
	#Check sign structure#
#	tm_exact,tp_exact = sign_structure(psi,Nsites,Nspec)
#
#	ft=open(res_path*"tau.txt","a")
#	println(ft,round(g,digits=4),"   ",tm_exact,"  ",tp_exact)
#	close(ft)
	
#	tm_sample,tp_sample = sample_sign(psi,Nsites,10000)
#	ft=open(res_path*"tau_sample.txt","a")
#	println(ft,round(g,digits=4),"   ",tm_sample,"  ",tp_sample)
#	close(ft)

#	Nsample=10000
#	real_part,imag_part = sample_sign_cplane(psi,Nsites,Nsample)
#	ft=open(res_path*"cplane_sample.txt","a")
#	for ii=1:Nsample
#		println(ft,round(g,digits=4),"   ",real_part[ii],"  ",imag_part[ii])
#	end
#	println(ft)
#	println(ft)
#	close(ft)
#	coeff =  sample_cplane(psi,Nsites,Nspec)
#	ft=open(res_path*"cplane.txt","a")
#	for ii=1:Nspec^3
#		println(ft,round(g,digits=4),"   ",real(coeff[ii]),"  ",imag(coeff[ii]))
#	end
#	println(ft)
#	println(ft)
#	close(ft)
	#################################################################


	mps_out=h5open(string(res_path*"psi0/"*"psi0_g",string(round(g,digits=3))),"w")
	write(mps_out,"MPS",psi)
	close(mps_out)

	#random_sample=sample(psi)

	global psi0 = psi

        f=open(res_path*"log","a")
        println(f," ")
        maxbond=maxlinkdim(psi)
        log_println(f,"Max. bond dimension: ",maxbond)
        println(f," ")
        @printf(f,"Final energy = %.8f \n",energy)
        println(f," ")
        close(f)

	if Nstates != 0

		energies = []
                append!(energies,energy)
                wavefunction = [psi for ii=1:Nstates+1]
                if ig == 0
                        global initial_states = [psi0 for ii=1:Nstates]
                end

		for istates=1:Nstates
			if ig == 0
	                	global initial_states[istates] = randomMPS(sites,Nbonds)
	        	else
	                	maxdim!(sweeps,maxlinkdim(initial_states[istates]),Nbonds) # gradually increase states kept
	        	end
	        	cutoff!(sweeps,SVD_error) # desired truncation error
	        	energy ,psi = dmrg(H,wavefunction[1:istates],initial_states[istates] ,sweeps, observer=obs, outputlevel=0)
	        	global initial_states[istates] = psi
	        	f2=open(res_path*"log","a")
	        	maxbond=maxlinkdim(psi)
	        	log_println(f,"Max. bond dimension: ",maxbond)
	        	println(f2)
	        	println(f2,"Final energy "*string(istates)*". excited state= "*string(round(energy,digits=12))*"\n")
        		println(f2)
        		close(f2)
        		wavefunction[istates+1]=psi
        		append!(energies,energy)
		end

		text_energy=" "
		text_ent_vN=" "
		text_ent_R=" "
		text_xcorr=" "
		text_ycorr=" "
		text_corr=" "
		text_mcorr=" "
		text_mux=" "
		text_muy=" "
		for istates=1:Nstates+1
			
			#Calculate von-Neumann entropy and Schmidt coefficients#
			SvN,Renyi,Svalues = vN_entropy(wavefunction[istates],mbond)
		
			write_output(res_path*"schmidt_values_"*string(istates)*".txt",g,Svalues)
			
			#Calculate dipole correlations#
			xcorr,ycorr = correlation(wavefunction[istates],Nsites,Nspec,evod,X,Y)
		
			#Write angular nearest-neighbour correlation to file#
			mcorr = ang_correlation(psi,Nsites,Nspec,evod,Lz)
			text_mcorr*=string((mcorr)/(Nsites-1)," ")

			#Calculate summed dipole moment and fluctuation#
			MuX,MuY = polarization(wavefunction[istates],Nsites,Nspec,evod,X,Y)
	
			#Write phi distribution#	
			Pphi = get_pphi(psi,Nsites,Nphi,Nspec,D)
			
			fd=open(res_path*"phi_distribution_"*string(istates)*".txt","a")
	                for ip=1:Nphi
	                        text=' '
	                        for b=1:Nsites
	                                text*=string(Pphi[b,ip],' ')
	                        end
	                        println(fd,round(g,digits=5),"   ",round(phi[ip],digits=4)," ",text)
	                end
	                println(fd)
	                println(fd)
	                close(fd)
		
			
			text_energy*=string(energies[istates]," ")
			text_ent_vN*=string(SvN," ")
			text_ent_R*=string(Renyi," ")
			text_xcorr*=string(xcorr," ")
			text_ycorr*=string(ycorr," ")
			text_corr*=string((xcorr+ycorr)/(Nsites-1)," ")
			text_mux*=string(MuX," ")
			text_muy*=string(MuY," ")

		end	

		write_text(res_path*"energy.txt",g,text_energy)
		write_text(res_path*"entropy.txt",g,text_ent)
		write_text(res_path*"xcorr.txt",g,text_xcorr)
		write_text(res_path*"ycorr.txt",g,text_ycorr)
		write_text(res_path*"corr.txt",g,text_corr)
		write_text(res_path*"mux.txt",g,text_mux)
		write_text(res_path*"muy.txt",g,text_muy)
		if evod == "all" || evod == "all_real"	
			write_text(res_path*"mcorr.txt",g,text_corr)
		end
	else
		#Write energy to file#
		f=open(res_path*"energy.txt","a")
		println(f,round(g,digits=4)," ",energy)
		close(f)
	
		#Write entropy and Schmidt values to file#
		SvN,Renyi,Svalues = vN_entropy(psi,mbond)
	
		write_output(res_path*"entropy_vN.txt",g,SvN)
		write_output(res_path*"entropy_Renyi.txt",g,Renyi)
		write_output(res_path*"schmidt_values.txt",g,Svalues)
		
		#Write nearest-neighbour correlation to file#
		xcorr,ycorr = correlation(psi,Nsites,Nspec,evod,X,Y)
		
		write_output(res_path*"xcorr.txt",g,xcorr)
		write_output(res_path*"ycorr.txt",g,ycorr)
		write_output(res_path*"corr.txt",g,(xcorr+ycorr)/(Nsites-1))

		#Write angular nearest-neighbour correlation to file#
		mcorr = ang_correlation(psi,Nsites,Nspec,evod,Lz)
		
		write_output(res_path*"mcorr.txt",g,(mcorr)/(Nsites-1))
		
		#Write polarization to file#
		MuX,MuY = polarization(psi,Nsites,Nspec,evod,X,Y)
	
		write_output(res_path*"mux.txt",g,MuX)
		write_output(res_path*"muy.txt",g,MuY)

		#Binder parameter#
		Mx2,My2,bindX,bindY = binder(psi,Nspec,X,Y)

		write_output(res_path*"binder_x.txt",g,[Mx2,bindX])
		write_output(res_path*"binder_y.txt",g,[My2,bindY])


		#Write phi distribution#	
		Pphi = get_pphi(psi,Nsites,Nphi,Nspec,D)

		fd=open(res_path*"phi_distribution.txt","a")
                for ip=1:Nphi
                        text=' '
                        for b=1:Nsites
                                text*=string(Pphi[b,ip],' ')
                        end
                        println(fd,round(g,digits=5),"   ",round(phi[ip],digits=4)," ",text)
                end
                println(fd)
                println(fd)
                close(fd)

	end

end
end
f=open(res_path*"log","a")
println(f," ")
log_println(f,"Calculation finished.")
close(f)

mv("log", res_path*"logDMRG")
cp("input.txt", res_path*"input.txt")