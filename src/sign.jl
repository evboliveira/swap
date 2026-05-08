module sign

using ITensors

export sign_structure,sample_sign,sample_sign_cplane,sample_cplane
##########################################################################################
function sign_structure(wf,Nsites,Nspec)

coeff=zeros(ComplexF64,(Nspec^3))
#coeff=zeros(Float64,(Nspec^6))

ii=1
for i1=1:Nspec
	s1=siteind(wf,1)
	proj1 = ITensor(s1)
	proj1[s1 => i1] = 1.0
	M1 = wf[1]*proj1
	for i2=1:Nspec
		s2=siteind(wf,2)
		proj2 = ITensor(s2)
		proj2[s2 => i2] = 1.0
		M2 = wf[2]*proj2
		for i3=1:Nspec
			s3=siteind(wf,3)
			proj3 = ITensor(s3)
			proj3[s3 => i3] = 1.0
			M3 = wf[3]*proj3
       
       			c=M1*M2*M3
       			coeff[ii] = scalar(c)
       			ii+=1
		end	
	end
end

tau_minus=0.0
tau_plus=0.0
for ii=1:Nspec^3
	c2 = real(coeff[ii])#*conj(coeff[ii]))
	if c2 < 0#coeff[ii] < 0
		tau_minus+=c2^2#coeff[ii]^2
	else
		tau_plus+=c2^2#coeff[ii]^2
	end		
end 

return tau_minus,tau_plus
end
##########################################################################################
function sample_sign(wf,Nsites,Nsample)

tminus = 0.0
tplus = 0.0
for jj=1:Nsample
	indices=sample(wf)

	s=siteind(wf,1)
	proj = ITensor(s)
	n=indices[1]
	proj[s => n] = 1.0
	psi = wf[1]*proj

	for ii=2:Nsites
		s=siteind(wf,ii)
		proj = ITensor(s)
		n=indices[ii]
		proj[s => n] = 1.0
		psi *= wf[ii]*proj
	end
	coeff=real(scalar(psi))
	if coeff < 0
		tminus+=1#coeff^2
	else
		tplus+=1#coeff^2
	end		
end
tau_minus = tminus/(tminus+tplus)
tau_plus = tplus/(tminus+tplus)

return tau_minus,tau_plus
end
##########################################################################################
function sample_sign_cplane(wf,Nsites,Nsample)

real_part=[]
imag_part=[]
for jj=1:Nsample
	indices=sample(wf)

	s=siteind(wf,1)
	proj = ITensor(s)
	n=indices[1]
	proj[s => n] = 1.0
	psi = wf[1]*proj

	for ii=2:Nsites
		s=siteind(wf,ii)
		proj = ITensor(s)
		n=indices[ii]
		proj[s => n] = 1.0
		psi *= wf[ii]*proj
	end
	coeff=real(scalar(psi))
	push!(real_part,real(coeff))
	push!(imag_part,imag(coeff))
end

return real_part,imag_part 
end
##########################################################################################
function sample_cplane(wf,Nsites,Nspec)

coeff=zeros(ComplexF64,(Nspec^3))

ii=1
for i1=1:Nspec
	s1=siteind(wf,1)
	proj1 = ITensor(s1)
	proj1[s1 => i1] = 1.0
	M1 = wf[1]*proj1
	for i2=1:Nspec
		s2=siteind(wf,2)
		proj2 = ITensor(s2)
		proj2[s2 => i2] = 1.0
		M2 = wf[2]*proj2
		for i3=1:Nspec
			s3=siteind(wf,3)
			proj3 = ITensor(s3)
			proj3[s3 => i3] = 1.0
			M3 = wf[3]*proj3
       
       			c=M1*M2*M3
       			coeff[ii] = scalar(c)
       			ii+=1
		end	
	end
end


return coeff 
end
##########################################################################################
end
