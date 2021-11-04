# -*- coding: utf-8 -*-
reflect <- function(incidentDir, normalVec){
  incidentDir - 2*c(crossprod(normalVec, incidentDir)) * normalVec
}

# n: number of segments; P0: initial point; v0: initial direction
trajectory <- function(n, P0, v0){
  out <- vector("list", n)
  L <- Line$new(P0, P0+v0)
  inters <- intersectionEllipseLine(ell, L)
  Q0 <- inters$I2
  out[[1]] <- Line$new(inters$I1, inters$I2, FALSE, FALSE)
  for(i in 2:n){
    theta <- atan2(Q0[2], Q0[1])
    t <- ell$theta2t(theta, degrees = FALSE)
    nrmlVec <- ell$normal(t)
    v <- reflect(Q0-P0, nrmlVec)
    inters <- intersectionEllipseLine(ell, Line$new(Q0, Q0+v))
    out[[i]] <- Line$new(inters$I1, inters$I2, FALSE, FALSE)
    P0 <- Q0
    Q0 <- if(isTRUE(all.equal(Q0, inters$I1))) inters$I2 else inters$I1
  }
  out
}

ell <- Ellipse$new(c(0,0), 6, 3, 0)

P0 <- ell$pointFromAngle(60)
v0 <- c(cos(pi+0.8), sin(pi+0.8))
traj <- trajectory(150, P0, v0)

opar <- par(mar = c(0,0,0,0))
plot(NULL, asp = 1, xlim = c(-7,7), ylim = c(-4,4),
     xlab = NA, ylab = NA, axes = FALSE)
draw(ell, border = "red", col = "springgreen", lwd = 3)
invisible(lapply(traj, draw))
par(opar)
```

Run the code below to see an animated trajectory:

```{r, eval=FALSE}
opar <- par(mar = c(0,0,0,0))
plot(NULL, asp = 1, xlim = c(-7,7), ylim = c(-4,4),
     xlab = NA, ylab = NA, axes = FALSE)
draw(ell, border = "red", col = "springgreen", lwd = 3)
for(i in 1:length(traj)){
  draw(traj[[i]])
  Sys.sleep(0.3)
}
