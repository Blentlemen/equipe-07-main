"""
Document qui donne les fonctions nécessaire pour le bon fonctionnement du PM de la question 3d
    - Place les dynodes sur le PM 12 dynodes précisément
    - Calcule le potentiel de chaque point du PM
    - Calcule le champ électrique du PM
    - Calcule le chemin de l'électron dans le PM avec les rebonds
    - Affiche le chemin sur le PM avec le potentiel et le champ électrique
"""

# Frédéric-Alexandre Caouette NI : 537275485
# Alexis Gauthier NI : 537289373
# Théo Mailhot NI : 537291728

import numpy as np


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


def champ_elec(res, parametres):
    """
    Fonction pour trouver le champ électrique du potentiel trouvé par la relaxation
        - Renvoie le champ électrique sous forme de matrice avec ses composantes (x, y)
        ainsi que la norme en tout point
    """

    dx = parametres[7]

    ey, ex = np.gradient(-res, dx, dx)
    e_norm = np.sqrt(ex**2 + ey**2)

    return [ex, ey, e_norm]


def eulerer_lechamp(x_p, y_p, ex, ey, dx):
    """
    Fonction pour trouver le champ en un point (x, y) précis
        - Renvoie les composantes spécifiques du champ électrique (x, y)
    """

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


def position_dynodes_bas(i, parametres):
    """
    Fonctions qui donnent les positions des dynodes en mm dans le PM
        - Renvoie les coordonnées (y_start:y_end, x_start:x_end) de chaque dynode
    """

    a, b, c, d, e, f = parametres[0:6]

    vert_start = -f*0.5 + b
    vert_end = vert_start + e
    horiz_start = a + i*(c + d)
    horiz_end = horiz_start + c
    pot = (2*i+1) * 100

    return [vert_start, vert_end, horiz_start, horiz_end, pot]

def position_dynodes_haut(i, parametres):
    """
    Fonctions qui donnent les positions des dynodes en mm dans le PM
        - Renvoie les coordonnées (y_start:y_end, x_start:x_end) de chaque dynode
    """

    a, b, c, d, e, f = parametres[0:6]

    vert_start = f*0.5 - b
    vert_end = vert_start - e
    horiz_start = a + (i+1)*c +d/2 + i*d - c/2
    horiz_end = horiz_start + c
    pot = (2*(i+1)) * 100

    return [vert_start, vert_end, horiz_start, horiz_end, pot]


def contact_dyn_bas(x_new, y_new, x_old, y_old, parametres):
    """
    Fonctions qui vérifient si l'électron entre en contact avec les dynodes
        - Renvoie un booléen pour savoir si il y a contact avec les dynode
          et renvoie la position de contact si True
    """

    f, n = parametres[5:7]
    lx = parametres[8]

    if -f/2 < y_new < 0 and 0 < x_new < lx:

        pente = (y_new - y_old) / (x_new - x_old)
        y_initiale = -(pente*x_old - y_old)

        for i in range (n//2 + n%2):

            dynodes_bas = position_dynodes_bas(i, parametres)

            if y_new <= dynodes_bas[1]:

                if x_new <= dynodes_bas[3] and x_old >= dynodes_bas[2]:

                    x_dyn = (dynodes_bas[1] - y_initiale)/pente
                    y_dyn = dynodes_bas[1]

                    return [True, x_dyn, y_dyn]

                if x_new >= dynodes_bas[2] and x_old <= dynodes_bas[2]:

                    if (dynodes_bas[0] - y_initiale)/pente >= dynodes_bas[2]:

                        if (dynodes_bas[1] - y_initiale)/pente <= dynodes_bas[2]:

                            x_dyn = dynodes_bas[2]
                            y_dyn = dynodes_bas[1]

                            return [True, x_dyn, y_dyn]

                        if (dynodes_bas[1] - y_initiale)/pente >= dynodes_bas[2]:

                            x_dyn = (dynodes_bas[1] - y_initiale)/pente
                            y_dyn = dynodes_bas[1]

                            return [True, x_dyn, y_dyn]

                if x_new >= dynodes_bas[3] and x_old <= dynodes_bas[3]:

                    if (dynodes_bas[1] - y_initiale)/pente <= dynodes_bas[3]:

                        x_dyn = (dynodes_bas[1] - y_initiale)/pente
                        y_dyn = dynodes_bas[1]

                        return  [True, x_dyn, y_dyn]

    return [False, x_new, y_new]

def contact_dyn_haut(x_new, y_new, x_old, y_old, parametres):
    """
    Fonctions qui vérifient si l'électron entre en contact avec les dynodes
        - Renvoie un booléen pour savoir si il y a contact avec les dynode
          et renvoie la position de contact si True
    """

    f, n = parametres[5:7]
    lx = parametres[8]

    if f/2 > y_new > 0 and 0 < x_new < lx:

        pente = (y_new - y_old) / (x_new - x_old)
        y_initiale = -(pente*x_old - y_old)

        for i in range (n//2):

            dynodes_haut = position_dynodes_haut(i, parametres)

            if y_new >= dynodes_haut[1]:

                if x_new <= dynodes_haut[3] and x_old >= dynodes_haut[2]:

                    x_dyn = (dynodes_haut[1] - y_initiale)/pente
                    y_dyn = dynodes_haut[1]

                    return [True, x_dyn, y_dyn]

                if x_new >= dynodes_haut[2] and x_old <= dynodes_haut[2]:

                    if (dynodes_haut[0] - y_initiale)/pente >= dynodes_haut[2]:

                        if (dynodes_haut[1] - y_initiale)/pente >= dynodes_haut[2]:

                            x_dyn = (dynodes_haut[1] - y_initiale)/pente
                            y_dyn = dynodes_haut[1]

                            return [True, x_dyn, y_dyn]

                        if (dynodes_haut[1] - y_initiale)/pente <= dynodes_haut[2]:

                            x_dyn = dynodes_haut[2]
                            y_dyn = dynodes_haut[1]

                            return [True, x_dyn, y_dyn]

                if x_new >= dynodes_haut[3] and x_old <= dynodes_haut[3]:

                    if (dynodes_haut[1] - y_initiale)/pente <= dynodes_haut[3]:

                        x_dyn = (dynodes_haut[1] - y_initiale)/pente
                        y_dyn = dynodes_haut[1]

                        return [True, x_dyn, y_dyn]

    return [False, x_new, y_new]


def position_el_rebond(x0, y0, vx0, vy0, ex, ey, dt, it_max, parametres):
    """
    Fonction qui donnent le chemin de l'électron dans le PM avec les rebonds sur les dynodes
        - Renvoie les coordonnées (x, y) de l'électron à chaque pas de temps dans un array
    """

    dx, lx, ly = parametres[7:10]

    q = -1.602*10e-19
    m = 9.109*10e-31

    x = [x0]
    y = [y0]
    vx = vx0
    vy = vy0

    ex_val, ey_val = eulerer_lechamp(x[0], y[0], ex, ey, dx)

    ax = (q*ex_val) / m
    ay = (q*ey_val) / m

    vx += ax * dt
    vy += ay * dt

    y_new = float(y0 + (vy * dt))
    x_new = float(x0 + (vx * dt))

    x.append(x_new)
    y.append(y_new)

    for _ in range(it_max):

        contact_bas = contact_dyn_bas(x_new, y_new, x[_], y[_], parametres)
        contact_haut = contact_dyn_haut(x_new, y_new, x[_], y[_], parametres)

        if contact_bas[0] is True:

            x_new = contact_bas[1]
            y_new = contact_bas[2]

            x.append(x_new)
            y.append(y_new)

            print(f"Rebond dynode BAS à [x,y] = [{x_new}, {y_new}]")
            print("La vitesse a été inversé et la position a AUGMENTÉ de 2mm")

            y_new += 2
            vy = -vy

            x.append(x_new)
            y.append(y_new)

        if contact_haut[0] is True:

            x_new = contact_haut[1]
            y_new = contact_haut[2]

            x.append(x_new)
            y.append(y_new)

            print(f"Rebond dynode HAUT à [x,y] = [{x_new}, {y_new}]")
            print("La vitesse a été inversé et la position a DIMINUÉ de 2mm")

            y_new -= 2
            vy = -vy

            x.append(x_new)
            y.append(y_new)

        ex_val, ey_val = eulerer_lechamp(x_new, y_new, ex, ey, dx)

        ax = (q*ex_val) / m
        ay = (q*ey_val) / m

        vx += ax * dt
        vy += ay * dt

        y_new = float(y_new + (vy * dt))
        x_new = float(x_new + (vx * dt))

        x.append(x_new)
        y.append(y_new)

        if x_new > lx and abs(y_new) < ly/2:
            print("L'électron à passé le photomultiplicateur au complet :)")
            break

        if not (0 <= x_new < lx and abs(y_new) < ly/2):
            print("L'électron a crissé son camp")
            break

    return np.array(x), np.array(y)
