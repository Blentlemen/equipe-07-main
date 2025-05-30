"""
Question 3b du travail numérique
"""

# Frédéric-Alexandre Caouette NI : 537275485
# Alexis Gauthier NI : 537289373
# Théo Mailhot NI : 537291728

import matplotlib.pyplot as plt
import question3b_fct

PARAMETRES = question3b_fct.initialiser_geo()

X0 = 0 # position initiale en mm
Y0 = 0 # position initiale en mm
VX0 = 0 # vitesse initiale en mm/s
VY0 = 0 # vitesse initiale en mm/s
DT = 1e-9 # pas de temps en s
IT_MAX = 10000 # nombre d'itérations maximum

X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

V, BLOQUER = question3b_fct.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
V, BLOQUER = question3b_fct.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

RES, DYNODES = question3b_fct.relaxation(V, BLOQUER, PARAMETRES, variation=1e-5, max_iter=10000)

EX, EY, E_NORM = question3b_fct.champ_elec(RES, PARAMETRES)[0:3]

traj_x, traj_y = question3b_fct.position_el(X0, Y0, VY0, VY0, EX, EY, DT, IT_MAX, PARAMETRES)

CP = plt.contourf(X_GRID, Y_GRID, RES, levels=100, cmap='plasma')
plt.contour(X_GRID, Y_GRID, BLOQUER, levels=[0.5], colors='black', linewidths=1)

plt.colorbar(CP, label="Potentiel (V)")
plt.plot(traj_x, traj_y, 'y-', label="Trajectoire")
plt.plot(X0, Y0, 'go', label="Départ (0,0)")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.title("Q3 b)")
plt.legend()
plt.axis("equal")
plt.savefig("Question3b)_trajectoire.png", dpi=300)
plt.show()

print("Position finale :", traj_x[-1], traj_y[-1])
print("Déplacement total :", traj_x[-1] - traj_x[0], traj_y[-1] - traj_y[0])
