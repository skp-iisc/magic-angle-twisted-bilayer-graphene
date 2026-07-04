import numpy as np
import matplotlib.pyplot as plt

th_deg = 0.5      # INPUT: twist in degrees

th = th_deg*np.pi/180
w = 110         # INPUT
d = 1.42          # angstrom
hv = 1.5*d*2970     # meV*angstrom
N = 5              # truncate range
valley = -1        # +1 for K, -1 for K'

KDens = 50     # INPUT: density of k points

kD = 8*np.pi*np.sin(th/2)/(3*d*np.sqrt(3))
b1m = kD*np.sqrt(3)*np.array([0.5, -np.sqrt(3)/2])
b2m = kD*np.sqrt(3)*np.array([0.5, np.sqrt(3)/2])
qb = kD*np.array([0, -1])
K1 = kD*np.array([-np.sqrt(3)/2, -0.5])
K2 = kD*np.array([-np.sqrt(3)/2, 0.5])

e1 = np.cos(2*np.pi/3) + complex(0,1)*valley*np.sin(2*np.pi/3)
e2 = np.cos(2*np.pi/3) - complex(0,1)*valley*np.sin(2*np.pi/3)

T_qb = w * np.array([[1, 1], [1, 1]], dtype=complex)
T_qtr = w * np.array([[e1, 1], [e2, e1]], dtype=complex)
T_qtl = w * np.array([[e2, 1], [e1, e2]], dtype=complex)
T_qb_dag = np.conjugate(T_qb.T)
T_qtr_dag = np.conjugate(T_qtr.T)
T_qtl_dag = np.conjugate(T_qtl.T)

# 2D Lattice
Lat = []
invL = np.zeros((2*N+1, 2*N+1), dtype=int)
count = 0
for i in range(-N, N+1):
    for j in range(-N, N+1):
        Lat.append([i, j])
        invL[i+N, j+N] = count
        count += 1
Lat = np.array(Lat)
# print(Lat, invL)
N_site = (2*N+1)**2

def h_k(kx, ky, ix, iy, Kval, valley):
    h_mat = np.zeros((2, 2), dtype=complex)
    if Kval == 'K1':
        ax = kx - valley*K1[0] + ix*b1m[0] + iy*b2m[0]
        ay = ky - valley*K1[1] + ix*b1m[1] + iy*b2m[1]
        qx = np.cos(th/2)*ax + np.sin(th/2)*ay
        qy = -np.sin(th/2)*ax + np.cos(th/2)*ay
    elif Kval == 'K2':
        ax = kx - valley*K2[0] + ix*b1m[0] + iy*b2m[0]
        ay = ky - valley*K2[1] + ix*b1m[1] + iy*b2m[1]
        qx = np.cos(th/2)*ax - np.sin(th/2)*ay
        qy = np.sin(th/2)*ax + np.cos(th/2)*ay
    h_mat[0, 1] = hv * (valley*qx - complex(0,1)*qy)
    h_mat[1, 0] = hv * (valley*qx + complex(0,1)*qy)
    return h_mat

def Ham_k_eigvals(kx, ky):
    Ham = np.zeros((4*N_site, 4*N_site), dtype=complex)
    for i in range(N_site):
        ix, iy = Lat[i]
        Ham[2*i, 2*i+1] = h_k(kx, ky, ix, iy, 'K1', valley)[0, 1]
        Ham[2*i+1, 2*i] = h_k(kx, ky, ix, iy, 'K1', valley)[1, 0]

        j = i + N_site
        Ham[2*j, 2*i] = T_qb_dag[0, 0]
        Ham[2*j, 2*i+1] = T_qb_dag[0, 1]
        Ham[2*j+1, 2*i] = T_qb_dag[1, 0]
        Ham[2*j+1, 2*i+1] = T_qb_dag[1, 1]

        if iy != valley*N:
            j = invL[ix+N, iy+valley+N] + N_site
            Ham[2*j, 2*i] = T_qtr_dag[0, 0]
            Ham[2*j, 2*i+1] = T_qtr_dag[0, 1]
            Ham[2*j+1, 2*i] = T_qtr_dag[1, 0]
            Ham[2*j+1, 2*i+1] = T_qtr_dag[1, 1]
        if ix != -valley*N:
            j = invL[ix-valley+N, iy+N] + N_site
            Ham[2*j, 2*i] = T_qtl_dag[0, 0]
            Ham[2*j, 2*i+1] = T_qtl_dag[0, 1]
            Ham[2*j+1, 2*i] = T_qtl_dag[1, 0]
            Ham[2*j+1, 2*i+1] = T_qtl_dag[1, 1]
    
    for i in range(N_site, 2*N_site):
        j = i - N_site
        ix, iy = Lat[j]
        Ham[2*i, 2*i+1] = h_k(kx, ky, ix, iy, 'K2', valley)[0, 1]
        Ham[2*i+1, 2*i] = h_k(kx, ky, ix, iy, 'K2', valley)[1, 0]

        Ham[2*j, 2*i] = T_qb[0, 0]
        Ham[2*j, 2*i+1] = T_qb[0, 1]
        Ham[2*j+1, 2*i] = T_qb[1, 0]
        Ham[2*j+1, 2*i+1] = T_qb[1, 1]
        if iy != -valley*N:
            j = invL[ix+N, iy-valley+N]
            Ham[2*j, 2*i] = T_qtr[0, 0]
            Ham[2*j, 2*i+1] = T_qtr[0, 1]
            Ham[2*j+1, 2*i] = T_qtr[1, 0]
            Ham[2*j+1, 2*i+1] = T_qtr[1, 1]
        if ix != valley*N:
            j = invL[ix+valley+N, iy+N]
            Ham[2*j, 2*i] = T_qtl[0, 0]
            Ham[2*j, 2*i+1] = T_qtl[0, 1]
            Ham[2*j+1, 2*i] = T_qtl[1, 0]
            Ham[2*j+1, 2*i+1] = T_qtl[1, 1]

    eig_vals, _ = np.linalg.eig(Ham)

    return np.sort(eig_vals)
