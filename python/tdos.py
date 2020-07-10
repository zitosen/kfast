from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

v = Vasprun('vasprun.xml')
tdos = v.tdos
plotter = DosPlotter(zero_at_efermi=True, stack=False, sigma=0.02)
plotter.add_dos("", tdos)
plotter.get_plot(xlim=[-3, 3], ylim=[0, 300]).savefig(fname="tdos.png",img_format="png")
