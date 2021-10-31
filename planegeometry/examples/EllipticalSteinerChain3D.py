# -*- coding: utf-8 -*-
import planegeometry.geometry as geom
import pyvista as pv
import numpy as np

c0 = geom.Circle((3,0), 3) # exterior circle
ell0 = geom.Ellipse((-4,0), 4, 2.5, 140)
f = geom.Affine.from_ellipse_to_ellipse(c0, ell0)

Tr = np.hstack(
    (
     np.vstack((f.A, np.zeros((2, 2), dtype=float))),
     np.column_stack(
         (
             np.array([0., 0., 1/3, 0.]),
             np.array([f.b[0], f.b[1], 0., 1.])
         )
     )
    )
)

# mesh = pv.Sphere(3, center=(3,0,0))
# mesh.transform(Tr)
# mesh.plot()

chain = geom.SteinerChain(c0, 3, -0.2, shift=50)
circles = chain["circles"]
ellipse_at_center = f.transform_ellipse(circles[-1])
cx, cy = ellipse_at_center.center
Center = (cx, cy, 0)
pltr = pv.Plotter()
for circle in circles:
    x, y = circle.center
    mesh = pv.Sphere(
        circle.radius, center=(x, y, 0), 
        phi_resolution=90, theta_resolution=90
    )
    mesh.transform(Tr)
    pltr.add_mesh(
        mesh, smooth_shading=True, color="firebrick", specular=10
    )
pltr.set_focus(Center)
pltr.set_position((0, 0, 15))
pltr.show()

