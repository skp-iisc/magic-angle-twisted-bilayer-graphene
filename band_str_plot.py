from band_str_calc import *
# CHANGE th_deg in band_str_calc.py to 1.05

A_to_B = np.arange(-1/2, 1/2, 1/KDens)
B_to_C = np.arange(-1, 0, 1/KDens)
C_to_D = np.arange(0, np.sqrt(3), 1/KDens)
D_to_A = np.arange(0, 1, 1/KDens)

K_len = len(A_to_B) + len(B_to_C) + len(C_to_D) + len(D_to_A)

Ek = []

for k in A_to_B:
    Ek.append(np.real(Ham_k_eigvals(np.sqrt(3)/2*kD, k*kD)))
for k in B_to_C:
    Ek.append(np.real(Ham_k_eigvals(0, k*kD)))
for k in C_to_D:
    Ek.append(np.real(Ham_k_eigvals(1/2*k*kD, -np.sqrt(3)/2*k*kD)))
for k in D_to_A:
    Ek.append(np.real(Ham_k_eigvals(-np.sqrt(3)/2*k*kD, -1/2*k*kD)))
Ek = np.array(Ek)

ylim = 80     # INPUT

plt.figure(figsize=(8,6))
for j in range(0, 4*N_site):
    plt.plot(np.arange(K_len), Ek[:,j], linestyle="-", linewidth=2)
plt.title(f"Moir$\\'e$ bands of twisted bilayer graphene, $\\theta={th_deg}^\\circ$",
        fontsize=14)
plt.xlim(0, K_len)
plt.ylim(-ylim, ylim)
plt.xticks([0, len(A_to_B), len(A_to_B)+len(B_to_C), 
            len(A_to_B)+len(B_to_C)+len(C_to_D), K_len], 
            ('A', 'B', 'C', 'D', 'A'), fontsize=12)
plt.vlines([0, len(A_to_B), len(A_to_B)+len(B_to_C), 
            len(A_to_B)+len(B_to_C)+len(C_to_D), K_len], 
            -ylim, ylim, colors='gray', 
            linestyles='dashed', linewidth=0.7)
plt.ylabel("Energy (meV)", fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('figs/fig3a_3.png', dpi=300)
plt.show()
