from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter

v = Vasprun('vasprun.xml')
cdos = v.complete_dos
element_dos = cdos.get_element_dos()
plotter = DosPlotter()
plotter.add_dos_dict(element_dos)
plotter.save_plot('pdos.png', img_format='png', xlim=[-3, 3], ylim=[0, 150])
