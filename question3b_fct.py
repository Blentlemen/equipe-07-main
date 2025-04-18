import numpy as np

"""
Fonction qui définit la géométrie du PM pour les question 1, 2, 3a) et 3b)
    - Renvoie tout les termes nécessaires pour les calculs des questions
"""
def initialiser_geo():

    #initialiser la géométrie
    a, b, c, d, e, f = 3, 2, 4, 2, 0.2, 6
    n = 4
    dx = 0.1 # step

    largeur = (2*a + (n+1)*(c/2) + (n-1)*(d/2))
    lx,ly = largeur, f
    nx, ny = int(lx/dx), int(ly/dx)

    x = np.linspace(0, lx, nx)
    y = np.linspace(-ly/2, ly/2, ny)
    x_grid, y_grid = np.meshgrid(x, y)

    #initialiser le potentiel à 0 partout
    v = np.zeros((ny,nx))
    bloquer = np.zeros((ny,nx), dtype=bool)

    return [a,b,c,d,e,f,n,dx,lx,ly,nx,ny,x,y,x_grid,y_grid,v,bloquer]

"""
Fonctions qui permettent de placer les dynodes sur le PM
    - Renvoie la matrice de potentiel avec les dynodes avec une
      matrice unique pour les dynodes qui empèche leur modification pendant la relaxation
"""
def placer_dynodes_bas(v, bloquer, parametres):

    a, b, c, d, e = parametres[0:5]
    n, dx = parametres[6:8]

    for i in range(n//2 + n%2):
        pot_dyn = (2*i+1) * 100
        vert_start = b
        vert_end = b + e
        horiz_start = a + i*(c+d)
        horiz_end = horiz_start + c

        ix_start = int(horiz_start/dx)
        ix_end = int(horiz_end/dx)
        iy_start = int(vert_start/dx)
        iy_end = int(vert_end/dx)

        v[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloquer[iy_start:iy_end, ix_start:ix_end] = True
    return v, bloquer

def placer_dynodes_haut(v, bloquer, parametres):

    a, b, c, d, e, f, n, dx = parametres[0:8]

    for i in range(n//2):

        pot_dyn = (2*(i+1)) * 100
        vert_start = f-b
        vert_end = vert_start + e
        horiz_start = a + (i+1)*c +d/2 + i*d - c/2
        horiz_end = horiz_start + c

        ix_start = int(horiz_start/dx)
        ix_end = int(horiz_end/dx)
        iy_start = int(vert_start/dx)
        iy_end = int(vert_end/dx)

        v[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloquer[iy_start:iy_end, ix_start:ix_end] = True

    return v, bloquer

"""
Fonction de relaxation pour trouver le potentiel en tout point
    - Renvoie le potentiel partout dans le PM
"""
def relaxation(v, bloquer, parametres, variation=1e-5, max_iter=1000000):

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

"""
Fonction pour trouver le champ électrique du potentiel trouvé par la relaxation
    - Renvoie le champ électrique sous forme de matrice avec ses composantes (x, y)
      ainsi que la norme en tout point
"""
def champ_elec(res, parametres):

    dx = parametres[7]

    ey, ex = np.gradient(-res, dx, dx)
    e_norm = np.sqrt(ex**2 + ey**2)

    return [ex, ey, e_norm]

"""
Fonction pour trouver le champ en un point (x, y) précis
    - Renvoie les composantes spécifiques du champ électrique (x, y)
"""
def eulerer_lechamp(x_p, y_p, ex, ey, dx):

    i = int(y_p / dx)
    j = int(x_p / dx)

    ny, nx = ex.shape

    ix = 0

    if i <= 0:

        ix = int(i + ny*0.5)

    if i > 0:

        ix = int(i + ny*0.5)

    if 0 <= ix <= ny and 0 <= j <= nx:

        return ex[ix, j], ey[ix, j]

    else:

        return 0, 0

"""
Fonction qui donne le déplacement d'un électron dans le PM
    - Renvoie les coordonnées (x, y) de l'électron à chaque pas de temps dans un array
"""
def position_el(x_ini, y_ini, vx_ini, vy_ini, ex, ey, step, max_it, parametres):

    q = -1.602e-19
    m = 9.109e-31

    x = [x_ini]
    y = [y_ini]
    vx = vx_ini
    vy = vy_ini

    dx = parametres[7]
    lx, ly = parametres[8:10]

    ex_val, ey_val = eulerer_lechamp(x_ini, y_ini, ex, ey, dx)

    ax = (q*ex_val) / m
    ay = (q*ey_val) / m

    vx += ax * step
    vy += ay * step

    x_new = float(x_ini + vx * step)
    y_new = float(y_ini + vy * step)

    x.append(x_new)
    y.append(y_new)

    for _ in range(max_it):

        ex_val, ey_val = eulerer_lechamp(x_new, y_new, ex, ey, dx)

        ax = (q*ex_val) / m
        ay = (q*ey_val) / m

        vx += ax * step
        vy += ay * step

        x_new = float(x_new + vx * step)
        y_new = float(y_new + vy * step)

        x.append(x_new)
        y.append(y_new)

        if not (0 <= x_new < lx and abs(y_new) <= ly/2):

            print("L'électron a crissé son camp")
            break

    print("Champ au point de départ :", eulerer_lechamp(x_ini, y_ini, ex, ey, dx))

    return np.array(x), np.array(y)
