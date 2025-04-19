"""
Question 1 du travail numérique
"""

# Frédéric-Alexandre Caouette NI : 537275485
# Alexis Gautier NI : 537289373
# Théo mailhot NI : 537291728

import matplotlib.pyplot as plt
import question1_fct


PARAMETRES = question1_fct.initialiser_geo()

X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

V, BLOQUER = question1_fct.placer_dynodes_haut(V, BLOQUER, PARAMETRES)

V, BLOQUER = question1_fct.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

RES, BLOQUER = question1_fct.relaxation(V, BLOQUER, PARAMETRES, variation=1e-5, max_iter=10000)
CP = plt.contourf(X_GRID, Y_GRID, RES, levels=100, cmap="plasma")

plt.contour(X_GRID, Y_GRID, BLOQUER, levels=[0.5], colors='black', linewidths=1)
plt.colorbar(CP, label="Potentiel (V)")
plt.title("Potentiel électrique dans le tube PM")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.axis('equal')
plt.tight_layout()
plt.savefig("Question1_potentiel.png", dpi=300)
plt.show()
