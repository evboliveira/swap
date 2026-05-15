

Nlist = [i for i=2:150]
listg = append!(
    [g for g = 0.0 : 0.1 : 0.4],
    [g for g = 0.41 : 0.01 : 0.60],
    [g for g = 0.7 : 0.1 : 3.0])
Ng = length(listg)

for Nsites=149:150
    psi0_path = "./results/psi0/"
    est_path = "./results/estimators/"
    if !isdir(psi0_path*"N$Nsites")
        mkdir(psi0_path*"N$Nsites")
    end
    if !isdir(est_path*"N$Nsites")
        mkdir(est_path*"N$Nsites")
    end
    res_path = "/home/evbdeoli/scratch/results/swap/N$Nsites/"

    for ig = 1:Ng
        g = listg[ig]
        f = string(res_path*"psi0/psi0_g",string(round(g,digits=3)))
        println(f)
        cp(f,string(psi0_path*"N$Nsites/psi0_g",string(round(g,digits=3))); force=true)
    end

    for file in readdir(res_path)
    # Check if the file has a .txt extension
        if endswith(file, ".txt")
            cp(joinpath(res_path, file), joinpath(est_path*"N$Nsites", file); force=true)
        end
    end
end
######################