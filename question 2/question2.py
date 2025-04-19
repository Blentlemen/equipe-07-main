"""
Question 2 du travail numérique
"""

# Frédéric-Alexandre Caouette NI : 537275485
# Alexis Gautier NI : 537289373
# Théo mailhot NI : 537291728

import matplotlib.pyplot as plt
import question2_fct


PARAMETRES = question2_fct.initialiser_geo()

X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

V, BLOQUER = question2_fct.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
V, BLOQUER = question2_fct.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

RES, BLOQUER = question2_fct.relaxation(V, BLOQUER, PARAMETRES, variation=1e-5, max_iter=10000)

EX, EY, E_NORM = question2_fct.champ_elec(RES, PARAMETRES)

plt.figure(figsize=(10, 5))
plt.title("Champ électrique dans le tube photomultiplicateur")
plt.contourf(X_GRID, Y_GRID, E_NORM, levels=100, cmap='plasma')
plt.colorbar(label="|E| (V/m)")

SAUT = 2

plt.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT],
            EX[::SAUT, ::SAUT], EY[::SAUT, ::SAUT],
            color='white', scale=10000)
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.axis('equal')
plt.tight_layout()
plt.savefig("Question2_champelec.png", dpi=300)
plt.show()
