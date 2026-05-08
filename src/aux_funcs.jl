module aux_funcs

export transform_realbasis, trans_matrix, create_file, write_output, write_output_error, write_text, log_println

###########     Module to store auxiliary functions for the main.jl file     ##########
#######################################################################################
function log_println(io_file, text...)
    println(stdout, text...)  # Print to terminal
    println(io_file, text...) # Print to file
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


########################### end of the module ###########################################
end