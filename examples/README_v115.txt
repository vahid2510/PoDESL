VTK EXPORT (called from Python after running dispatch):
from podesl.solvers.vtk_writer import write_vtk_quad2d
res = dispatch(problem, env)
write_vtk_quad2d("ps_result.vtk", res["nodes"], res["elems"],
                 point_data={
                    "U": res["U"].reshape(-1,2),
                    "sigma_x": res["sigma"][:,0],
                    "sigma_y": res["sigma"][:,1],
                    "tau_xy":  res["sigma"][:,2],
                    "ep":      res["ep"]
                 })
