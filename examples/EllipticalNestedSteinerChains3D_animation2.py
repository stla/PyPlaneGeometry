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

def add_ellipse(pltr, circle, color):
    x, y = circle.center
    mesh = pv.Sphere(
        circle.radius, center=(x, y, 0), 
        phi_resolution=90, theta_resolution=90
    )
    mesh.transform(Tr)
    pltr.add_mesh(
        mesh, smooth_shading=True, color=color, specular=10
    )

palettes = ["viridis", "magma", "inferno", "plasma", "cividis"]

shift_ = np.linspace(0, 5, 101)[:100]
for i, shift in enumerate(shift_):
    circles = geom.SteinerChain(c0, 5, -0.15, shift)["circles"]
    if i == 0:
        ellipse_at_center = f.transform_ellipse(circles[-1])
        cx, cy = ellipse_at_center.center
        Center = (cx, cy, 0)
    # Steiner chain for each circle, except the small one (it is too small)
    chains = [
        geom.SteinerChain(c, 3, -0.4, -shift*3/5)["circles"] for c in circles[:5]
    ]
    pltr = pv.Plotter(window_size=[512, 512], off_screen=True)
    pltr.set_background("#363940")
    for j, chain in enumerate(chains):
        colors = sb.color_palette(palette=palettes[j], n_colors=4)
        for k, circle in enumerate(chain):
            add_ellipse(pltr, circle, colors[k])
    pltr.set_focus(Center)
    pltr.set_position((0, -15, 15))
    pltr.camera.zoom(1.5)
    pngname = "zpic_%03d.png" % i
    pltr.show(screenshot=pngname)

os.system(
    "magick convert -dispose previous -loop 0 -delay 8 zpic_*.png EllipticalNestedSteinerChains3D_2.gif"    
) 
