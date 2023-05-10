# -*- coding: utf-8 -*-
import os
import pyvista as pv
import planegeometry.geometry as geom
import numpy as np
import seaborn as sb


c0 = geom.Circle((7,0), 3) # exterior circle
ell0 = geom.Ellipse((0,0), 4, 2.5, 140)
f = geom.Affine.from_ellipse_to_ellipse(c0, ell0)

# 4x4 transformation matrix
Tr = np.hstack(
    (
     np.vstack((f.A, np.zeros((2, 2), dtype=float))),
     np.column_stack(
         (
             np.array([0., 0., 0.75, 0.]),
             np.array([f.b[0], f.b[1], 0., 1.])
         )
     )
    )
)

def add_ellipse(pltr, circle, Center):
    x, y = circle.center
    mesh = pv.Sphere(
        circle.radius, center=(x, y, 0), 
        phi_resolution=90, theta_resolution=90
    )
    mesh.transform(Tr)
    Centers = np.repeat(Center.reshape((1,3)), len(mesh.points), axis=0)
    mesh.point_data["distance"] = np.linalg.norm(mesh.points-Centers, axis=1)
    pltr.add_mesh(
        mesh, smooth_shading=True, cmap="rocket", specular=10, 
        show_scalar_bar=False
    )


shift_ = np.linspace(0, 5, 101)[:100]
for i, shift in enumerate(shift_):
    circles = geom.SteinerChain(c0, 5, -0.15, shift)["circles"]
    if i == 0:
        ellipse_at_center = f.transform_ellipse(circles[-1])
        cx, cy = ellipse_at_center.center
        Center = np.array([cx, cy, 0])
    # Steiner chain for each circle, except the small one (it is too small)
    chains = [
        geom.SteinerChain(c, 3, -0.2, -shift*3/5)["circles"] for c in circles[:5]
    ]
    pltr = pv.Plotter(window_size=[512, 512], off_screen=True)
    pltr.set_background("#363940")
    for chain in chains:
        for circle in chain:
            add_ellipse(pltr, circle, Center)
    pltr.set_focus(Center)
    pltr.set_position((0, -15, 15))
    pltr.camera.zoom(1.4)
    pngname = "zpic_%03d.png" % i
    pltr.show(screenshot=pngname)

os.system(
    "magick convert -dispose previous -loop 0 -delay 8 zpic_*.png EllipticalNestedSteinerChains3D_3.gif"    
) 
