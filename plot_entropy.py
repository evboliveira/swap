import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
from matplotlib import gridspec
from matplotlib import ticker
from matplotlib.patches import Patch
from matplotlib import cm
from matplotlib.ticker import LinearLocator
# matplotlib.use("TkAgg")

#### deffining colors
blue = (40/255, 95/255, 185/255)
red = (195/255, 58/255, 50/255)
green = (0/255, 158/255, 115/255)
orange = (230/255, 159/255, 0/255)

# Nlist = [5,10,50,100,150]
Nlist = [i for i in range(20,150,10)]
nN = len(Nlist)
glist = np.round(np.append(np.linspace(0,0.4,5), np.append(np.linspace(0.41,0.59,19),np.linspace(0.6,3,25))), decimals=3)
Ng = len(glist)

for n in Nlist:
    ## Data extraction DMRG ##
    res_path = "/home/evbdeoli/scratch/results/swap/20may2026/N%i/"%n
    # res_path = "/home/evbdeoli/scratch/results/swap/20may2026/N%i/"%Nlist[iN]
    # swap = np.loadtxt(res_path+"entropy_swap.txt")[0:,1]
    # swap_err = np.loadtxt(res_path+"entropy_swap.txt")[0:,2]
    NMswap = np.loadtxt(res_path+"NMentropy_swap.txt")[0:,1]
    NMswap_err = np.loadtxt(res_path+"NMentropy_swap.txt")[0:,2]
    plt.plot(glist,NMswap, label='N=%i'%n)
    plt.xlabel("g")
    plt.ylabel("Ent. Entropy")
    plt.legend()
plt.savefig('entropy_vs_g.png')
plt.close()

Nlist = [i for i in range(5,150)]
nN = len(Nlist)
list = np.zeros(nN)
for g in [0.1,0.5,1.0, 2.0, 3.0]:
    print(np.where(glist==g))
    ig = np.where(glist==g)[0][0]
    for iN in range(nN):
        ## Data extraction DMRG ##
        res_path = "/home/evbdeoli/scratch/results/swap/20may2026/N%i/"%Nlist[iN]
        # res_path = "/home/evbdeoli/scratch/results/swap/20may2026/N%i/"%Nlist[iN]
        list[iN] = np.loadtxt(res_path+"NMentropy_swap.txt")[ig,1]
        
    plt.plot(Nlist, list, label='g=%.2f'%g)
    plt.xlabel("N")
    plt.ylabel("Ent. Entropy")
    plt.legend()
plt.savefig('entropy_vs_N.png')
plt.close()


######################
## 3D plot #####

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


Nlist = np.array([i for i in range(4,50,1)])
nN = len(Nlist)
glist = np.round(np.append(np.linspace(0,0.4,5), np.append(np.linspace(0.41,0.59,19),np.linspace(0.6,3,25))), decimals=3)
Ng = len(glist)

data = np.zeros((Ng,nN))
for iN in range(nN):
    ## Data extraction DMRG ##
    res_path = "/home/evbdeoli/scratch/results/swap/20may2026/N%i/"%Nlist[iN]
    # res_path = "/home/evbdeoli/scratch/results/swap/20may2026/N%i/"%Nlist[iN]
    data[:,iN] = np.loadtxt(res_path+"NMentropy_swap.txt")[0:,1]

# def Z_func(g,N):
#     return data[np.where(glist==g)[0][0], np.where(Nlist==N)[0][0]]

X, Y = np.meshgrid(Nlist,glist)
# print(Y)
# Z = Z_func(X,Y)
surf = ax.plot_surface(X, Y, data, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)


# # Customize the z axis.
ax.set_zlim(0, 0.5)
ax.zaxis.set_major_locator(LinearLocator(10))
# # A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.set_xlabel("N")
ax.set_ylabel("g")
ax.set_zlabel("Ent. Entropy")
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig("grid.png")
plt.close()


## Phase Diagram ##
plt.imshow(data,
           interpolation='gaussian',
           cmap=cm.coolwarm,
           extent=(Nlist[0], Nlist[-1], glist[-1], glist[0]),
           aspect='auto')


plt.savefig("phase.png")