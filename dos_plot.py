from band_str_calc import *

ks = np.arange(0, 1, 1/KDens)
K_len = len(ks)**2
Ek = []
for m in ks:
    for n in ks:
        kx = (m-n)*np.sqrt(3)/2
        ky = (m+n)*1.5 -1
        Ek.append(np.real(Ham_k_eigvals(kx*kD, ky*kD)))
Ek = np.array(Ek)
E_vals = Ek.flatten()

xlim = 50    #INPUT
plt.figure(figsize=(8,6))
plt.hist(E_vals, bins=5000, color='blue', alpha=0.7)
plt.title(f"Density of States ($\\theta={th_deg}^\\circ$)",
        fontsize=16)
plt.xlabel("Energy (meV)", fontsize=14)
plt.ylabel("Density of States", fontsize=14)
plt.xlim(-xlim, xlim)
# E_plt = []
# for E in E_vals:
#     if abs(E) < xlim:
#         E_plt.append(E)
# E_plt = np.array(E_plt)
# ylim = max(E_plt)*1.1
# plt.ylim(0, ylim)
plt.yticks([])
plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig('figs/fig3b_3.png', dpi=300)
plt.show()
