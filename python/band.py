import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter

v = BSVasprun("vasprun.xml")
bs = v.get_band_structure(kpoints_filename="KPOINTS",line_mode=True)
band_fig = BSPlotter(bs)
band_fig.get_plot(vbm_cbm_marker=True,ylim=(-1.5,3.5))
plt.savefig('band.png', img_format='png')
