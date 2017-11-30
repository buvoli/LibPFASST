import glob
import numpy as np
from os import remove, mkdir
from nose import with_setup
from pf.pfasst import PFASST, Experiment, Params

# defined relative to root of project
exe = 'tests/toda/main.exe'
base_dir = 'tests/toda/output'
try:
    mkdir(base_dir)
except:
    pass


params = Params(nodes=[3], magnus=[2], tasks=1, periodic=True, \
                nsteps=128, tfinal=1.0, particles=11, iterations=15, \
                nb=False, base_dir=base_dir, solutions=False)

toda = PFASST(exe, params) 
results = toda.run()[0]
final_solution = results.loc[len(results)-1, 'solution']
periodic_ref = final_solution

params = Params(nodes=[3], magnus=[2], tasks=1, periodic=False, \
                nsteps=128, tfinal=1.0, particles=11, iterations=15, \
                nb=False, base_dir=base_dir, solutions=False)

toda = PFASST(exe, params) 
results = toda.run()[0]
final_solution = results.loc[len(results)-1, 'solution']
nonperiodic_ref = final_solution

def cleanup():
    for fname in glob.iglob(base_dir+'/*pkl'):
        remove(fname)

def test_generator():
    for periodic in [True, False]:
        for nodes in [[3], [2]]:
            for magnus in [[2], [1]]:
                if nodes[0] == 2 and magnus[0] == 2: continue
                for tasks in [1]:
                    for tol in [0.0, 1e-15]:
                        yield toda, nodes, magnus, tasks, periodic, tol

@with_setup(teardown=cleanup)
def toda(nodes, magnus, tasks, periodic, tol):
    params = Params(nodes=nodes, magnus=magnus, periodic=periodic, 
                    tasks=tasks, tolerance=tol,
                    nsteps=128, tfinal=1.0, particles=11, iterations=30, 
                    nb=False, base_dir=base_dir, solutions=False)

    toda = PFASST(exe, params) 
    results = toda.run()[0]
    final_solution = results.loc[len(results)-1, 'solution']

    if periodic:
        ref = periodic_ref
    else:
        ref = nonperiodic_ref

    print 'max number of iterations: {}'.format(results.iter.max())
    print 'mean residual: {}'.format(results.residual.mean())
    print 'max difference: {}'.format(np.amax(final_solution - ref))
    assert np.amax(final_solution - ref) < 1e-5