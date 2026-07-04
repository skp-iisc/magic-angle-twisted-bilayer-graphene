# Magic-angle twisted bilayer graphene (MA-TBG)

Solution of BM Hamiltonian.

## Band structure plots
Choose `th_deg` and take `KDens = 100` in  the `band_struc_calc.py` file.
- `th_deg = 5`, ylim(-1600, 1600).
- `th_deg = 1.05`, ylim(-300, 300).
- `th_deg = 0.5`, ylim(-80, 80).

## Density of states plots
Choose `th_deg` and take `KDens = 50` in  the `band_struc_calc.py` file.
- `th_deg = 5`, xlim(-1000, 1000).
- `th_deg = 1.05`, xlim(-200, 200).
- `th_deg = 0.5`, xlim(-50, 50).

----------

## Other details

Lat: [[-5, -5], [-5, -4], ..., [-5, 5], [-4, -5], ..., [5, 5]] for N=5.

invL: $7\times7$ 2d array for N=3, with values from 0 to 48. $11\times11$ 2d array for N=5, with values from 0 to 120.

