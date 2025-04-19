"""
Question 3d du travail numérique
"""

import numpy as np
import matplotlib.pyplot as plt
import question3d_fct

A, B, C, D, E, F = 4, 3, 6, 4, 0.2, 10
N = 12
DX = 0.1

LARGEUR = (2*A + (N+1)*(C/2) + (N-1)*(D/2))
LX,LY = LARGEUR, F
NX, NY = int(LX/DX), int(LY/DX)

X = np.linspace(0, LX, NX)
Y = np.linspace(-LY/2, LY/2, NY)
X_GRID, Y_GRID = np.meshgrid(X, Y)

V = np.zeros((NY, NX))
BLOQUER = np.zeros((NY, NX), dtype=bool)

PARAMETRES = A,B,C,D,E,F,N,DX,LX,LY,NX,NY,X,Y,X_GRID,Y_GRID,V,BLOQUER

V, BLOQUER = question3d_fct.placer_dynodes_bas(V, BLOQUER, PARAMETRES)
V, BLOQUER = question3d_fct.placer_dynodes_haut(V, BLOQUER, PARAMETRES)

RES, BLOQUER = question3d_fct.relaxation(V, BLOQUER, PARAMETRES, variation=1e-5, max_iter=10000)

EX, EY, E_NORM = question3d_fct.champ_elec(RES, PARAMETRES)[0:5]

X0 = 0
Y0 = 0
VX0 = 0
VY0 = 0
DT = 0.00000001
IT_MAX = 50000

traj_x, traj_y = question3d_fct.position_el_rebond(X0, Y0, VX0, VY0, EX, EY, DT, IT_MAX, PARAMETRES)


fig, axe = plt.subplots(figsize=(10, 6))

CP = axe.contourf(X_GRID, Y_GRID, E_NORM, levels=100, cmap='plasma')

plt.colorbar(CP, label="|E| (V/m)")

SAUT = 2
"""
axe.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT], EX[::SAUT, ::SAUT],
               EY[::SAUT, ::SAUT], color='white', scale=6000)

axe.set_xlim(0, LX)
axe.set_ylim(-LY/2, LY/2)
axe.set_xlabel("x (mm)")
axe.set_ylabel("y (mm)")
axe.set_title("Animation de la trajectoire de l'électron")
axe.legend()
axe.axis('equal')
"""
#-------------------------------affichage de la trajectoire x(t)-------------------------------#

plt.contourf(X_GRID, Y_GRID, E_NORM, levels=100, cmap='plasma')
plt.contour(X_GRID, Y_GRID, BLOQUER, levels=[0.5], colors='black', linewidths=1)

SAUT = 2

plt.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT],
        EX[::SAUT, ::SAUT], EY[::SAUT, ::SAUT],
        color='white', scale=6000)
plt.plot(traj_x, traj_y, 'y-', label="Trajectoire")
plt.plot(X0, Y0, 'go', label="Départ (x0,y0)")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.title("Trajectoire")
plt.legend()
plt.axis("equal")
plt.savefig("Question3c)_trajectoire.png", dpi=300)
plt.show()

print("Position finale :", traj_x[-1], traj_y[-1])
print("Déplacement total :", traj_x[-1] - traj_x[0], traj_y[-1] - traj_y[0])
print("Champ au point de départ :", question3d_fct.eulerer_lechamp(X0, X0, EX, EY, DX))
