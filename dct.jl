using FFTW
v=[0,1,2,-3,4,5,6,4,-8,5,7]
println(v)
u=dct(v)*sqrt(size(v)[1])
#println(u)
w=idct(u)/sqrt(size(v)[1])
for i=1:size(w)[1]
    w[i]=round(Int,w[i])
    w[i]=convert(Int,w[i])
end
println(w)
