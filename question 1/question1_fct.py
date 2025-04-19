"""
Document qui donne les fonctions nécessaires pour la question 1
    - Place les dynodes sur le PM
    - Calcule le potentiel en tout point du PM
    - Affiche le potentiel sur le PM
"""

# Frédéric-Alexandre Caouette NI : 537275485
# Alexis Gautier NI : 537289373
# Théo mailhot NI : 537291728

import numpy as np

def initialiser_geo():
    """
    Fonction qui définit la géométrie du PM pour les question 1, 2, 3a) et 3b)
        - Renvoie tout les termes nécessaires pour les calculs des questions
    """

    #initialiser la géométrie
    a, b, c, d, e, f = 3, 2, 4, 2, 0.2, 6
    n = 4
    dx = 0.1 # step

    largeur = (2*a + (n+1)*(c/2) + (n-1)*(d/2))
    lx,ly = largeur, f
    nx, ny = int(lx*(1/dx)), int(ly*(1/dx))

    x = np.linspace(0, lx, nx)
    y = np.linspace(-ly/2, ly/2, ny)
    x_grid, y_grid = np.meshgrid(x, y)

    #initialiser le potentiel à 0 partout
    v = np.zeros((ny,nx))
    bloquer = np.zeros((ny,nx), dtype=bool)

    return [a,b,c,d,e,f,n,dx,lx,ly,nx,ny,x,y,x_grid,y_grid,v,bloquer]


def placer_dynodes_bas(v, bloquer, parametres):
    """
    Fonctions qui permettent de placer les dynodes sur le PM
        - Renvoie la matrice de potentiel avec les dynodes avec une
        matrice unique pour les dynodes qui empèche leur modification pendant la relaxation
    """

    a, b, c, d, e = parametres[0:5]
    n, dx = parametres[6:8]

    for i in range (n//2 + n%2):
        a, b, c, d, e = parametres[0:5]
        dx = parametres[7]

        pot_dyn = (2*i+1) * 100
        vert_start = b
        vert_end = b + e
        horiz_start = a + i*(c+d)
        horiz_end = horiz_start + c

        ix_start = int(horiz_start*(1/dx))
        ix_end = int(horiz_end*(1/dx))
        iy_start = int(vert_start*(1/dx))
        iy_end = int(vert_end*(1/dx))

        v[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloquer[iy_start:iy_end, ix_start:ix_end] = True
    return v, bloquer

def placer_dynodes_haut(v, bloquer, parametres):
    """
    Fonctions qui permettent de placer les dynodes sur le PM
        - Renvoie la matrice de potentiel avec les dynodes avec une
        matrice unique pour les dynodes qui empèche leur modification pendant la relaxation
    """

    a, b, c, d, e, f, n, dx = parametres[0:8]

    for i in range (n//2):
        a, b, c, d, e, f = parametres[0:6]
        dx = parametres[7]

        pot_dyn = (2*(i+1)) * 100
        vert_start = f - b - e
        vert_end = vert_start + e
        horiz_start = a + (i + 1 - 0.5)*c + (i + 0.5)*d
        horiz_end = horiz_start + c

        ix_start = int(horiz_start*(1/dx))
        ix_end = int(horiz_end*(1/dx))
        iy_start = int(vert_start*(1/dx))
        iy_end = int(vert_end*(1/dx))

        v[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloquer[iy_start:iy_end, ix_start:ix_end] = True
    return v, bloquer


def relaxation(v, bloquer, parametres, variation=1e-5, max_iter=1000000):
    """
    Fonction de relaxation pour trouver le potentiel en tout point
        - Renvoie le potentiel partout dans le PM
    """

    nx, ny = parametres[10:12]

    for iteration in range(max_iter):

        v_old = v.copy()

        for i in range(1, ny-1):

            for j in range(1, nx-1):

                if not bloquer[i, j]:

                    v[i, j] = 0.25 * (v_old[i+1, j] + v_old[i-1, j] +
                                     v_old[i, j+1] + v_old[i, j-1])

        chang_pot_max = np.max(np.abs(v - v_old))

        if chang_pot_max <= variation:

            print(f"Convergence atteinte après {iteration} itérations "\
                  "avec une tolérance de variation de {variation}")
            print("Nous avons atteint un potentiel qui ne varie presque plus " \
            "et le problème est considéré résolu")

            break

        if iteration == max_iter - 1:

            print("Attention !!!!!!!!!! ")
            print("Le maximum d'itérations a été atteint sans stabilisation, " \
            "donc le programme a été arrêté")

    return v, bloquer
