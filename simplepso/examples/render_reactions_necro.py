"""Reaction graph for Figure 5C from the PySB publication"""

from __future__ import print_function
import pysb.tools.render_reactions
from pysb.kappa import contact_map, set_kappa_path, influence_map
from pysb.tools.render_reactions import run
from necroptosismodule import model
import pygraphviz as pyg

# print out the graphviz representation of the contact map
# print(pysb.tools.render_reactions.run(model))


set_kappa_path('/Users/geenaildefonso/Projects/KaSim')
# x = contact_map(model)
# x.draw('contact_map.pdf', format='pdf', prog='dot')

x = run(model)
g = pyg.AGraph(x)
g.draw('render_reactions_necro_vt_n.pdf', format='pdf', prog='dot')