import os
import matplotlib.pyplot as plt
import planegeometry.geometry as geom
import seaborn as sb
import numpy as np


# function to draw a circle
def draw_circle(axes, C, color, fill=True):
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=fill, facecolor=color, 
            edgecolor="black", linewidth=2
        )
    )

# function to construct the "children" ####
def children(inversions, circles1):
    m = len(inversions)
    n = len(circles1) 
    circles2 = [] 
    for i in range(n):
        k = circles1[i][1]
        for j in range(m):
            if j != k:
                circle = inversions[j].invert_circle(circles1[i][0])
                circles2.append((circle, j))
    return circles2

#
c0 = geom.Circle((0,0), 3) # the exterior circle
n = 3
depth = 5
colors = sb.color_palette(palette="bright", n_colors=depth)
shift_ = np.linspace(0, n, 101)[:100]

for k, shift in enumerate(shift_): 
    # starting circles ####
    circles0 = geom.SteinerChain(c0, n, phi = 0.25, shift = shift)["circles"]
    
    # construct the inversions ####
    inversions = [None]*(n+1)
    for i in range(n):
        inversions[i] = geom.Inversion.from_fixing_three_circles(
            c0, circles0[i], circles0[(i+1) % n]
        )
    inversions[n] = geom.Inversion.from_swapping_two_circles(
        c0, circles0[n]
    )
    
    # first generation of children
    circles1 = []
    for i in range(n):
        ip1 = (i+1) % n
        for j in range(n+1): #(j in (1L:(n+1L))[-c(i,ip1)]){
            if j != i and j != ip1:
                circle = inversions[i].invert_circle(circles0[j])
                circles1.append((circle, i))
        
    # construct children ####
    allCircles = [None]*depth
    allCircles[0] = circles0
    allCircles[1] = circles1
    for i in range(depth)[2:]:
        allCircles[i] = children(inversions, allCircles[i-1])
    for i in range(depth)[1:]:
        allCircles[i] = [c[0] for c in allCircles[i]]
    
    # plot ####
    figure, axes = plt.subplots(figsize=(10, 10))
    axes.set_aspect(1)
    draw_circle(axes, c0, "black", False)
    for i in range(depth):
        for circ in allCircles[i]:
            draw_circle(axes, circ, colors[i])
    plt.title("Apollonian gasket", fontdict = {"fontsize": 40})
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.axis("off")
    pngname = "zpic_%03d.png" % k
    plt.savefig(pngname, dpi=96)

os.system(
    "mogrify -resize 512x zpic_*.png"    
) 
os.system(
    "magick convert -dispose previous -loop 0 -delay 8 zpic_*.png ApollonianGasket.gif"    
) 