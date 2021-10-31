# -*- coding: utf-8 -*-
import os
import planegeometry.geometry as geom
import pyvista as pv
import numpy as np

c0 = geom.Circle((3,0), 3) # exterior circle
ell0 = geom.Ellipse((-4,0), 4, 2.5, 140)
f = geom.Affine.from_ellipse_to_ellipse(c0, ell0)

# 4x4 transformation matrix
Tr = np.hstack(
    (
     np.vstack((f.A, np.zeros((2, 2), dtype=float))),
     np.column_stack(
         (
             np.array([0., 0., 0.5, 0.]),
             np.array([f.b[0], f.b[1], 0., 1.])
         )
     )
    )
)

shift_ = np.linspace(0, 3, 101)[:100]
for i, shift in enumerate(shift_):
    chain = geom.SteinerChain(c0, 3, -0.2, shift=shift)
    circles = chain["circles"]
    if i == 0:
        ellipse_at_center = f.transform_ellipse(circles[-1])
        cx, cy = ellipse_at_center.center
        Center = (cx, cy, 0)
    pltr = pv.Plotter(window_size=[512, 512], off_screen=True)
    pltr.set_background("#363940")
    for circle in circles:
        x, y = circle.center
        mesh = pv.Sphere(
            circle.radius, center=(x, y, 0), 
            phi_resolution=90, theta_resolution=90
        )
        mesh.transform(Tr)
        pltr.add_mesh(
            mesh, smooth_shading=True, color="crimson", specular=10
        )
    pltr.set_focus(Center)
    pltr.set_position((0, -15, 15))
    pltr.camera.zoom(1.5)
    pngname = "zpic_%03d.png" % i
    pltr.show(screenshot=pngname)

os.system(
    "magick convert -dispose previous -loop 0 -delay 8 zpic_*.png EllipticalSteinerChain3D.gif"    
) 