ITensors.space(::SiteType"PlaRotor") = 11

#######################################################################
function ITensors.op!(Op::ITensor,::OpName"T",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = T[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"V6",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = V6[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"X",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = X[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"X2",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = X2[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"X3",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = X3[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"X4",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = X4[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"Y",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = Y[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"Y2",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = Y2[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"Y3",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = Y3[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"Y4",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = Y4[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"Ycomp",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = 1.0im*Y[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"Lz",::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = Lz[j,i]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D1" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,1]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D2" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,2]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D3" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,3]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D4" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,4]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D5" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,5]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D6" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,6]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D7" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,7]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D8" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,8]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D9" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,9]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D10" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,10]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D11" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,11]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D12" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,12]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D13" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,13]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D14" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,14]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D15" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,15]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D16" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,16]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D17" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,17]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D18" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,18]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D19" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,19]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D20" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,20]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D21" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,21]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D22" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,22]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D23" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,23]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D24" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,24]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D25" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,25]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D26" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,26]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D27" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,27]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D28" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,28]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D29" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,29]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D30" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,30]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D31" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,31]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D32" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,32]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D33" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,33]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D34" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,34]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D35" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,35]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D36" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,36]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D37" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,37]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D38" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,38]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D39" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,39]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D40" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,40]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D41" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,41]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D42" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,42]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D43" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,43]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D44" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,44]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D45" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,45]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D46" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,46]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D47" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,47]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D48" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,48]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D49" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,49]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D50" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,50]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D51" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,51]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D52" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,52]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D53" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,53]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D54" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,54]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D55" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,55]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D56" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,56]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D57" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,57]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D58" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,58]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D59" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,59]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D60" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,60]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D61" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,61]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D62" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,62]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D63" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,63]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D64" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,64]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D65" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,65]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D66" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,66]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D67" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,67]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D68" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,68]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D69" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,69]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D70" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,70]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D71" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,71]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D72" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,72]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D73" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,73]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D74" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,74]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D75" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,75]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D76" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,76]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D77" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,77]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D78" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,78]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D79" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,79]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D80" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,80]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D81" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,81]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D82" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,82]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D83" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,83]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D84" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,84]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D85" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,85]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D86" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,86]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D87" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,87]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D88" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,88]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D89" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,89]
end
end

end
#######################################################################

function ITensors.op!(Op::ITensor,::OpName"D90" ,::SiteType"PlaRotor" ,s::Index)

complex!(Op)
for i=1:Nspec
for j=1:Nspec
     Op[s'=>j,s=>i] = D[j,i,90]
end
end

end
#######################################################################

