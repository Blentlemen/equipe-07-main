import matplotlib.pyplot as plt
import fonctions


PARAMETRES = fonctions.initialiser_geo()

X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

V, BLOQUER = fonctions.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
V, BLOQUER = fonctions.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

RES, BLOQUER = fonctions.relaxation(V, BLOQUER, PARAMETRES, variation=1e-5, max_iter=10000)

EX, EY, E_NORM = fonctions.champ_elec(RES, PARAMETRES)

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
