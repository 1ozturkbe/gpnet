from __future__ import print_function
from __future__ import absolute_import
from gpkit.constraints.tight import Tight
from gpnet.draw_network import draw_KT_network
from gpnet.testing.layout import hanoi, hanoi_from_data
import numpy as np

from robust.robust import RobustModel

if __name__ == '__main__':
    # Changing Tightness requirement to within 0.1%
    Tight.reltol = 10.**-3
    m = hanoi(friction='HW')
    warm_start = {m["D"]: (1.016 - 0.3048)*np.random.rand(len(m.topology_list)) + 0.3048}
    sol = m.localsolve(verbosity=2, reltol=1e-2, iteration_limit=150, x0=warm_start)
    print(sol['cost'])

    draw_KT_network(sol, m.coordinates, m.topology_list)

    # Comparing against robust solution
    rm = RobustModel(m, 'elliptical', nominalsolve=sol, gamma=1)
    rmsol = rm.robustsolve()
    draw_KT_network(rmsol, m.coordinates, m.topology_list)
