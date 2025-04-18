import question3a_fct

PARAMETRES = question3a_fct.initialiser_geo()

X0 = 0
Y0 = 0
VX0 = 0
VY0 = 0
DT = 1e-9
IT_MAX = 10000

X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

V, BLOQUER = question3a_fct.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
V, BLOQUER = question3a_fct.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

RES, BLOQUER = question3a_fct.relaxation(V, BLOQUER, PARAMETRES, variation=1e-5, max_iter=10000)

EX, EY = question3a_fct.champ_elec(RES, PARAMETRES)[0:2]

traj_x, traj_y = question3a_fct.position_el(X0, Y0, VY0, VY0, EX, EY, DT, IT_MAX, PARAMETRES)