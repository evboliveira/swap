module create_operators

export operators
#########################################################################################
function operators(Nspec,Nphi,evod)

f=open("operators.jl","w")

#Define the Hilbert space#
println(f,"ITensors.space(::SiteType\"PlaRotor\") = ", Nspec)
println(f)

#Define string building blocks#
stringop="complex!(Op)"
string1="function ITensors.op!(Op::ITensor,::OpName"
string2=",::SiteType\"PlaRotor\" ,s::Index)"
string3="for i=1:Nspec"
string4="for j=1:Nspec"
string5="       Op[s'=>j,s=>i] = D[j,i,"
string6="end"
string7="#######################################################################"


#Write kinetic energy operator#
println(f,string7)
println(f,string1,"\"T\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = T[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#Write V6  operator#
println(f,string1,"\"V6\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = V6[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

###################
#Write X operators#
###################

#X#
println(f,string1,"\"X\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = X[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#X2#
println(f,string1,"\"X2\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = X2[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#X3#
println(f,string1,"\"X3\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = X3[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#X4#
println(f,string1,"\"X4\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = X4[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

###################
#Write Y operators#
###################

#Y#
println(f,string1,"\"Y\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = Y[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#Y2#
println(f,string1,"\"Y2\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = Y2[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#Y3#
println(f,string1,"\"Y3\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = Y3[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)

#Y4#
println(f,string1,"\"Y4\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = Y4[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)


println(f,string1,"\"Ycomp\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = 1.0im*Y[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)


#Write Lz operator#
println(f,string1,"\"Lz\"",string2)
println(f)
if evod == "all"
	println(f,stringop)
end
println(f,string3)
println(f,string4)
println(f,"     Op[s'=>j,s=>i] = Lz[j,i]")
println(f,string6)
println(f,string6)
println(f)
println(f,string6)
println(f,string7)
println(f)


#Write phi-distribution operator#
for ip=1:Nphi
	name = string("\"D",ip,"\" ")
	println(f,string1,name,string2)
	println(f)
	if evod == "all"
		println(f,stringop)
	end
	println(f,string3)
	println(f,string4)
	println(f,"     Op[s'=>j,s=>i] = D[j,i,",ip,"]")
	println(f,string6)
	println(f,string6)
	println(f)
	println(f,string6)
	println(f,string7)
	println(f)
end
close(f)

end
#########################################################################################
end
