import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
from scipy.special import ellipe, ellipk

from matplotlib import gridspec
from matplotlib import ticker

from matplotlib.patches import Patch
# matplotlib.use("TkAgg")

#### deffining colors
blue = (40/255, 95/255, 185/255)
red = (195/255, 58/255, 50/255)
green = (0/255, 158/255, 115/255)
orange = (230/255, 159/255, 0/255)

Nlist = [i for i in range(5,151)]
nN = len(Nlist)
glist = np.round(np.append(np.linspace(0,0.4,5), np.append(np.linspace(0.41,0.59,19),np.linspace(0.6,3,25))), decimals=3)
Ng = len(glist)

for n in Nlist:
    ## Data extraction DMRG ##
    res_path = "./results/N%i/"%n
    swap = np.loadtxt(res_path+"entropy_swap.txt")[0:,1]
    swap_err = np.loadtxt(res_path+"entropy_swap.txt")[0:,2]
    NMswap = np.loadtxt(res_path+"NMentropy_swap.txt")[0:,1]
    NMswap_err = np.loadtxt(res_path+"NMentropy_swap.txt")[0:,2]

# ##Plotting results ####
# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'''
#     \usepackage{nicefrac}
#     \usepackage{siunitx}
#     \usepackage{physics} 
#     \usepackage[mathcal]{eucal}
#     \DeclareFontFamily{OT1}{mathc}{}
#     \DeclareFontShape{OT1}{mathc}{m}{it}{<-> mathc10}{}
#     \DeclareMathAlphabet{\mathcal}{OT1}{mathc}{m}{it}
#     '''
# )
# plt.rc('font', family='serif',serif='Times')

# revtex_col = 3.39
# w = revtex_col
# h = 2.5
# figsize = (w,h)
# FontSize = 10
# cols = 2
# rows = 2

# tickpad=0.9
# LabelPad=2.5 #y axis


# fig, [[ax1,ax3],[ax2,ax4]] = plt.subplots(rows, cols, figsize=figsize, sharex=True)
# fig.subplots_adjust(hspace=0.0, wspace=0.35, left=0.115, right=0.997, top=0.999, bottom=0.16)

# ## ENERGY
# # fig.subplots_adjust(hspace=0.4,wspace=0.0)
# ax1.text(
#     0.45, 0.95, r"\bf(a)",
#     transform=ax1.transAxes,
#     fontsize=FontSize,
#     fontweight='bold',
#     va='top'
# )
# ax1.plot(G[:glimidx], Edmrg[:glimidx], zorder=1, lw=1.2, ls='-', color = 'black',label=r'$\varepsilon_{{\scriptscriptstyle \text{DMRG}}}$')

# ax1.plot(G[:idx_gc], Epert[:idx_gc], zorder=3, ls='-', lw=1.2, color = red,label=r'$\varepsilon_{\text{dis}}$')
# ax1.plot(G[idx_gc:idx_gc+dgc], Epert[idx_gc:idx_gc+dgc], zorder=2, ls='--', lw = 1, color = red,label=None)

# ax1.plot(G[0:idx_gc], Esaa[0:idx_gc], zorder=2, ls='--', lw = 1, color = blue,label=None)
# ax1.plot(G[idx_gc:glimidx], Esaa[idx_gc:glimidx], zorder=3, ls='-', lw=1.2, color = blue,label=r'$\varepsilon_{\text{ord}} - \frac{1}{8}$')

# ax1.vlines([gc], 0, 1, color='gray', linestyle='dashed', lw=1., label=None,zorder=1,
#     transform=ax1.get_xaxis_transform())
# # Plot parameters
# ax1.legend(fontsize=0.8*FontSize,loc='lower left', handlelength=1, labelspacing=0.15 )
# ax1.set_xlabel(r'$g$', fontsize = FontSize)
# ax1.set_ylabel(r'$\varepsilon(g)$', fontsize = FontSize, labelpad=LabelPad)#, rotation=0, labelpad=15)
# # ax1.minorticks_on()
# # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
# ax1.set_ylim(-2.4,0.4)
# ax1.set_xticks([0, gc, 1, 2])
# ax1.set_xticklabels([0, r'$\sim g_c ~~~~~~~$', 1, 2])
# ax1.xaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize)
# ax1.yaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize, pad=tickpad)
# plt.setp(ax1.spines.values(), linewidth=1)  


list = np.zeros(nN)
for g in [0.1,0.5,1.0, 2.0, 3.0]:
    print(np.where(glist==g))
    ig = np.where(glist==g)[0][0]
    for iN in range(nN):
        ## Data extraction DMRG ##
        # res_path = "./results/N%i/"%Nlist[iN]
        res_path = "/home/evbdeoli/scratch/results/swap/N%i/"%Nlist[iN]
        list[iN] = np.loadtxt(res_path+"NMentropy_swap.txt")[ig,1]
        
    plt.plot(Nlist, list, label='g=%.2f'%g)
    plt.xlabel("N")
    plt.ylabel("Ent. Entropy")
    plt.legend()
plt.savefig('test.png')


# res_file = "results//"

# #### Data extraction DMRG ####
# G = np.loadtxt(res_file+"energy.txt", delimiter='\t', dtype=float, encoding=None)[1:,0]
# idx_gc = np.where(G==gc)[0][0]
# dgc = 40

# Kdmrg = np.loadtxt(res_file+"KE.txt", delimiter='\t', dtype=float, encoding=None)[1:,1]/N
# Edmrg = np.loadtxt(res_file+'energy.txt', delimiter='\t', dtype=float, encoding=None)[1:,1]/N
# # V_dmrg = E_dmrg - K_dmrg
# L2dmrg = np.loadtxt(res_file+'Ltotal2.txt', delimiter='\t', dtype=float, encoding=None)[1:,1]/N
# # Mcorr_dmrg = np.loadtxt(res_file+'mcorr.txt', delimiter=' ', dtype=float, encoding=None)[:,1]
# binderx_dmrg = np.loadtxt(res_file+'binder_x.txt', delimiter='\t', dtype=float, encoding=None)[1:,2]

# gtest = np.loadtxt(res_file+'mux_copy.txt', delimiter='\t', dtype=float, encoding=None)[::-1,0]
# Xpol_dmrg_test = np.loadtxt(res_file+'mux_copy.txt', delimiter='\t', dtype=float, encoding=None)[::-1,1]/N
# # print(gtest)
# # print(Xpol_dmrg_test)

# phi_corr_dmrg = (np.loadtxt(res_file+'xcorr.txt', delimiter='\t', dtype=float, encoding=None)[1:,1]/N
#                  +
#                  np.loadtxt(res_file+'ycorr.txt', delimiter='\t', dtype=float, encoding=None)[1:,1]/N)

# len_dmrg = G.shape[0]
# plt_range_dmrg = range(1,len_dmrg-1)


# #### Small Angle Approximation  + phi4 ####
# Esaa = -2*G + (2*np.sqrt(3)*ellipe(2/3)/np.pi)*np.sqrt(G) -1/8
# Epert = -(5/8)*G**2

# Ediff_saa = np.absolute((Edmrg-Esaa))/Edmrg
# Ediff_pert = np.absolute((Edmrg-Epert))/Edmrg

# L2saa = 0.5*np.sqrt(3*G) -1/8 #- 0.1*G**(-1/2)
# L2pert = (9/8)*G**2
# L2diff_saa = np.absolute((L2dmrg-L2saa))/L2dmrg
# L2diff_pert = np.absolute((L2dmrg-L2pert))/L2dmrg

# cos_t = 1
# Msaa = cos_t*(1- (ellipk(2/3)/(2*np.pi*np.sqrt(3)))/np.sqrt(G))
# Mpert = np.array([0 for i in G])
# Mdiff_per = np.absolute((Mpert-Xpol_dmrg_test)/Xpol_dmrg_test)
# Mdiff_saa = np.absolute((Msaa-Xpol_dmrg_test)/Xpol_dmrg_test)

# binder_saa = 1 - (1/3)*(1-ellipk(2/3)/(np.pi*np.sqrt(3*G)))**(-1)

# Xcorr_saa = 1 + (np.sqrt(3)*(ellipe(2/3)-ellipk(2/3))/np.pi)/np.sqrt(G)
# Xcorr_per = (1/4)*G
# Xcoor_diff_per = np.absolute((Xcorr_per-phi_corr_dmrg)/phi_corr_dmrg)
# Xcoor_diff_saa = np.absolute((Xcorr_saa-phi_corr_dmrg)/phi_corr_dmrg)


# glimidx=200
# ##-----------------------------------------------------


# ##Plotting results ####
# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'''
#     \usepackage{nicefrac}
#     \usepackage{siunitx}
#     \usepackage{physics} 
#     \usepackage[mathcal]{eucal}
#     \DeclareFontFamily{OT1}{mathc}{}
#     \DeclareFontShape{OT1}{mathc}{m}{it}{<-> mathc10}{}
#     \DeclareMathAlphabet{\mathcal}{OT1}{mathc}{m}{it}
#     '''
# )
# plt.rc('font', family='serif',serif='Times')

# revtex_col = 3.39
# w = revtex_col
# h = 2.5
# figsize = (w,h)
# FontSize = 10
# cols = 2
# rows = 2

# tickpad=0.9
# LabelPad=2.5 #y axis


# fig, [[ax1,ax3],[ax2,ax4]] = plt.subplots(rows, cols, figsize=figsize, sharex=True)
# fig.subplots_adjust(hspace=0.0, wspace=0.35, left=0.115, right=0.997, top=0.999, bottom=0.16)

# ## ENERGY
# # fig.subplots_adjust(hspace=0.4,wspace=0.0)
# ax1.text(
#     0.45, 0.95, r"\bf(a)",
#     transform=ax1.transAxes,
#     fontsize=FontSize,
#     fontweight='bold',
#     va='top'
# )
# ax1.plot(G[:glimidx], Edmrg[:glimidx], zorder=1, lw=1.2, ls='-', color = 'black',label=r'$\varepsilon_{{\scriptscriptstyle \text{DMRG}}}$')

# ax1.plot(G[:idx_gc], Epert[:idx_gc], zorder=3, ls='-', lw=1.2, color = red,label=r'$\varepsilon_{\text{dis}}$')
# ax1.plot(G[idx_gc:idx_gc+dgc], Epert[idx_gc:idx_gc+dgc], zorder=2, ls='--', lw = 1, color = red,label=None)

# ax1.plot(G[0:idx_gc], Esaa[0:idx_gc], zorder=2, ls='--', lw = 1, color = blue,label=None)
# ax1.plot(G[idx_gc:glimidx], Esaa[idx_gc:glimidx], zorder=3, ls='-', lw=1.2, color = blue,label=r'$\varepsilon_{\text{ord}} - \frac{1}{8}$')

# ax1.vlines([gc], 0, 1, color='gray', linestyle='dashed', lw=1., label=None,zorder=1,
#     transform=ax1.get_xaxis_transform())
# # Plot parameters
# ax1.legend(fontsize=0.8*FontSize,loc='lower left', handlelength=1, labelspacing=0.15 )
# ax1.set_xlabel(r'$g$', fontsize = FontSize)
# ax1.set_ylabel(r'$\varepsilon(g)$', fontsize = FontSize, labelpad=LabelPad)#, rotation=0, labelpad=15)
# # ax1.minorticks_on()
# # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
# ax1.set_ylim(-2.4,0.4)
# ax1.set_xticks([0, gc, 1, 2])
# ax1.set_xticklabels([0, r'$\sim g_c ~~~~~~~$', 1, 2])
# ax1.xaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize)
# ax1.yaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize, pad=tickpad)
# plt.setp(ax1.spines.values(), linewidth=1)

# #### L2
# ax2.text(
#     0.45, 0.95, r"\bf(b)",
#     transform=ax2.transAxes,
#     fontsize=FontSize,
#     fontweight='bold',
#     va='top'
# )
# # Second plot: top-right (position 0,1 in a 2x2 grid)
# ax2.plot(G[:glimidx], L2dmrg[:glimidx], zorder=2, ls='-',lw=1.2, color = 'black',label=r'$\ell^2_{{\scriptscriptstyle \text{DMRG}}}$')
# ax2.plot(G[:idx_gc], L2pert[:idx_gc], zorder=3, ls='-', lw=1.2, color = red,label=r'$\ell^2_{\text{dis}}$')
# ax2.plot(G[idx_gc:idx_gc+dgc], L2pert[idx_gc:idx_gc+dgc], zorder=2, ls='--', lw = 1, color = red,label=None)
# ax2.plot(G[1:idx_gc], L2saa[1:idx_gc], zorder=2, ls='--', lw = 1, color = blue,label=None)
# ax2.plot(G[idx_gc:glimidx], L2saa[idx_gc:glimidx], zorder=3, ls='-', lw=1.2, color = blue,label=r'$\ell^2_{\text{ord}}-\frac{1}{8}$')
# ax2.vlines([gc], 0, 1, color='gray', linestyle='dashed', lw=1., label=None,zorder=1,transform=ax2.get_xaxis_transform())
# # Plot parameters
# ax2.legend(fontsize=0.8*FontSize,loc='lower right', handlelength=1, labelspacing=0.35 )
# ax2.set_xlabel(r'$g$', fontsize = FontSize)
# ax2.set_ylabel(r'$\ell^2(g)$', fontsize = FontSize, labelpad=LabelPad)#, rotation=0, labelpad=15)
# # ax2.minorticks_on()
# # ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
# ax2.set_ylim(-0.1,1.2)
# ax2.set_xticks([0, gc, 1, 2])
# ax2.set_xticklabels([0, r'$\sim g_c ~~~~~~~$', 1, 2])
# ax2.xaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize)
# ax2.yaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize, pad=tickpad)
# ax2.tick_params(axis='y', which='minor', width=0.5, length=1)
# ax2.tick_params(axis='x', which='minor', width=0.5, length=1)
# for i, label in enumerate(ax2.get_xticklabels()):
#     if i == 1:  # Change size of 'Two' (index 2)
#         label.set_fontsize(FontSize) # Optional: change color to highlight
#     else:
#         label.set_fontsize(0.8*FontSize)
# ax2.xaxis.set_tick_params(width=1, length=3)
# ax2.xaxis.get_major_ticks()[1].tick1line.set_markersize(5)
# ax2.xaxis.get_major_ticks()[1].tick1line.set_markeredgewidth(1.2)
# plt.setp(ax2.spines.values(), linewidth=1)


# #### Mag
# ax3.text(
#     0.45, 0.95, r"\bf(c)",
#     transform=ax3.transAxes,
#     fontsize=FontSize,
#     fontweight='bold',
#     va='top'
# )
# # Second plot: top-right (position 0,1 in a 2x2 grid)
# ax3.plot(gtest[:glimidx], Xpol_dmrg_test[:glimidx], zorder=2, ls='-',lw=1.2, color = 'black',label=r'$\mathcal{m}_{{\scriptscriptstyle \text{DMRG}}}$')
# ax3.plot(G[:idx_gc], Mpert[:idx_gc], zorder=3, ls='-', lw=1.2, color = red,label=r'$\mathcal{m}_{\text{dis}}$')
# ax3.plot(G[idx_gc:idx_gc+dgc], Mpert[idx_gc:idx_gc+dgc], zorder=2, ls='--', lw = 1, color = red,label=None)
# ax3.plot(G[10:idx_gc], Msaa[10:idx_gc], zorder=2, ls='--', lw = 1, color = blue,label=None)
# ax3.plot(G[idx_gc:glimidx], Msaa[idx_gc:glimidx], zorder=3, ls='-', lw=1.2, color = blue,label=r'$\mathcal{m}_{\text{ord}}$')
# ax3.vlines([gc], 0, 1, color='gray', linestyle='dashed', lw=1., label=None,zorder=1,transform=ax3.get_xaxis_transform())
# # Plot parameters
# ax3.legend(fontsize=0.8*FontSize,loc='lower right', handlelength=1, labelspacing=0.35 )
# ax3.set_xlabel(r'$g$', fontsize = FontSize)
# ax3.set_ylabel(r'$\mathcal{m}(g)$', fontsize = FontSize, labelpad=LabelPad)#, rotation=0, labelpad=15)
# # ax3.set_xlim(0,2)
# ax3.set_ylim(-0.1,1.1)
# ax3.set_yticks([0, 0.5, 1])
# ax3.set_yticklabels([r'$0.0$', r'$0.5$', r'$1.0$' ])
# # ax3.minorticks_on()
# # ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
# ax3.set_xticks([0, gc, 1, 2])
# ax3.set_xticklabels([0, r'$\bf \sim g_c ~~~~~~~$', 1, 2])
# ax3.xaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize)
# ax3.yaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize, pad=tickpad)
# ax3.tick_params(axis='y', which='minor', width=0.5, length=1)
# ax3.tick_params(axis='x', which='minor', width=0.5, length=1)
# plt.setp(ax3.spines.values(), linewidth=1)


# #### Corr
# ax4.text(
#     0.45, 0.95, r"\bf(d)",
#     transform=ax4.transAxes,
#     fontsize=FontSize,
#     fontweight='bold',
#     va='top'
# )
# # Second plot: top-right (position 0,1 in a 2x2 grid)
# ax4.plot(G[:glimidx], phi_corr_dmrg[:glimidx], zorder=2, ls='-',lw=1.2, color = 'black',label=r'$\mathcal{c}_{{\scriptscriptstyle \text{DMRG}}}$')
# ax4.plot(G[:idx_gc], Xcorr_per[:idx_gc], zorder=3, ls='-', lw=1.2, color = red,label=r'$\mathcal{c}_{\text{dis}}$')
# ax4.plot(G[idx_gc:idx_gc+dgc], Xcorr_per[idx_gc:idx_gc+dgc], zorder=2, ls='--', lw = 1, color = red,label=None)
# ax4.plot(G[30:idx_gc], Xcorr_saa[30:idx_gc], zorder=2, ls='--', lw = 1, color = blue,label=None)
# ax4.plot(G[idx_gc:glimidx], Xcorr_saa[idx_gc:glimidx], zorder=3, ls='-', lw=1.2, color = blue,label=r'$\mathcal{c}_{\text{ord}}$')
# ax4.vlines([gc], 0, 1, color='gray', linestyle='dashed', lw=1., label=None,zorder=1,transform=ax4.get_xaxis_transform())
# # Plot parameters
# ax4.legend(fontsize=0.8*FontSize,loc='lower right', handlelength=1, labelspacing=0.35 )
# ax4.set_xlabel(r'$g$', fontsize = FontSize)
# ax4.set_ylabel(r'$\mathcal{c}(g)$', fontsize = FontSize, labelpad=LabelPad)#, rotation=0, labelpad=15)
# ax4.set_ylim(-0.1,0.8)
# ax4.set_yticks([0,0.2, 0.4, 0.6])
# ax4.set_yticklabels([r'$0.0$', r'$0.2$',r'$0.4$', r'$0.6$' ])
# # ax4.minorticks_on()
# # ax4.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
# ax4.set_xticks([0, gc, 1, 2])
# ax4.set_xticklabels([0, r'$\bf \sim g_c$', 1, 2])
# for i, label in enumerate(ax4.get_xticklabels()):
#     if i == 1:  # Change size of 'Two' (index 2)
#         label.set_fontsize(FontSize) # Optional: change color to highlight
#     else:
#         label.set_fontsize(0.8*FontSize)
# ax4.xaxis.set_tick_params(width=1, length=3)
# ax4.xaxis.get_major_ticks()[1].tick1line.set_markersize(5)
# ax4.xaxis.get_major_ticks()[1].tick1line.set_markeredgewidth(1.2)
# ax4.yaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize, pad=tickpad)
# plt.setp(ax4.spines.values(), linewidth=1)# plt.tight_layout()
# plt.savefig("N=150_grid.png",format = 'png', dpi=600)



# revtex_col = 3.39
# w = revtex_col
# h = 1.4
# figsize = (w,h)
# FontSize = 10
# cols = 1
# rows = 1

# fig, ax5 = plt.subplots(rows, cols, figsize=figsize, sharex=True)
# fig.subplots_adjust(hspace=0.0, left=0.11, right=0.99, top=0.95, bottom=0.18)

# glimidx =290

# line1, = ax5.plot(G[:idx_gc], Ediff_pert[:idx_gc], zorder=2, ls='-', lw=1.2, color = blue,label=r'$\frac{\left|\varepsilon_{\text{dis/ord}}-\varepsilon_{{\scriptscriptstyle \text{DMRG}}}\right|}{\varepsilon_{\scriptscriptstyle \text{DMRG}}}$')
# ax5.plot(G[idx_gc:glimidx], Ediff_saa[idx_gc:glimidx], zorder=2, ls='-', lw=1.2, color = blue,label=None)

# line2, = ax5.plot(G[:idx_gc], L2diff_pert[:idx_gc], zorder=2, ls='-', lw=1.2, color = red,label=r'$\frac{\left|\ell^2_{\text{dis/ord}}-\ell^2_{{\scriptscriptstyle \text{DMRG}}}\right|}{\ell_{\scriptscriptstyle \text{DMRG}}}$')
# ax5.plot(G[idx_gc:glimidx], L2diff_saa[idx_gc:glimidx], zorder=2, ls='-', lw=1.2, color = red,label=None)

# line3, = ax5.plot(G[:idx_gc], Mdiff_per[:idx_gc], zorder=2, ls='-', lw=1.2, color = orange,label=r'$\frac{\left|\mathcal{m}_{\text{dis/ord}}-\mathcal{m}_{{\scriptscriptstyle \text{DMRG}}}\right|}{\mathcal{m}_{\scriptscriptstyle \text{DMRG}}}$')
# ax5.plot(G[idx_gc:glimidx], Mdiff_saa[idx_gc:glimidx], zorder=2, ls='-', lw=1.2, color = orange,label=None)

# line4, = ax5.plot(G[:idx_gc], Xcoor_diff_per[:idx_gc], zorder=2, ls='-', lw=1.2, color = green,label=r'$\frac{\left|\mathcal{c}_{\text{dis/ord}}-\mathcal{c}_{{\scriptscriptstyle \text{DMRG}}}\right|}{\mathcal{c}_{\scriptscriptstyle \text{DMRG}}}$')
# ax5.plot(G[idx_gc:glimidx], Xcoor_diff_saa[idx_gc:glimidx], zorder=2, ls='-', lw=1.2, color = green,label=None)

# ax5.hlines(0, 0, 1, color='black', linestyle='-', lw=1.2, label=None,zorder=1,
#     transform=ax5.get_yaxis_transform())
# ax5.vlines(gc, 0, 1, color='gray', linestyle='dashed', lw=1., label=None,zorder=1,
#     transform=ax5.get_xaxis_transform())
# # Plot parameters
# ax5.set_xscale('log')
# ax5.set_ylim(-0.25,0.75)
# ax5.legend(fontsize=0.8*FontSize,loc='best', handlelength=1, labelspacing=0.5 )
# ax5.set_xlabel(r'$g$', fontsize = FontSize)
# # axs.set_ylabel(r'Rel. Diff.', fontsize = FontSize)#, rotation=0, labelpad=15)
# # ax5.minorticks_on()
# ax5.set_xticks([1e-2, 1e-1, gc, 1e0, 1e1])
# ax5.set_xticklabels([r'$10^{-2}$',r'$10^{-1}$', r'$\bf\sim\!g_c$', r'$10^{0}$', r'$10^{1}$'])
# ax5.yaxis.set_tick_params(width=1, length=3, labelsize=0.8*FontSize, pad=tickpad)
# for i, label in enumerate(ax5.get_xticklabels()):
#     if i == 2:  # Change size of 'Two' (index 2)
#         label.set_fontsize(FontSize) # Optional: change color to highlight
#     else:
#         label.set_fontsize(0.8*FontSize)
# ax5.xaxis.set_tick_params(width=1, length=3)
# ax5.xaxis.get_major_ticks()[2].tick1line.set_markersize(5)
# ax5.xaxis.get_major_ticks()[2].tick1line.set_markeredgewidth(1.3)
# # Create first legend
# first_legend = ax5.legend(handles=[line1,line2], loc='upper left',fontsize = 0.95*FontSize, handlelength=1, labelspacing=0.5 )
# ax5.add_artist(first_legend)

# # Create second legend
# second_legend = ax5.legend(handles=[line3,line4], loc='upper right',fontsize = 0.95*FontSize, handlelength=1, labelspacing=0.5 )
# ax5.add_artist(second_legend)

# plt.setp(ax5.spines.values(), linewidth=1)

# # plt.tight_layout()
# plt.savefig("N=150_rel_diff.png",format = 'png', dpi=600)
