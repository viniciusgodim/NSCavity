library(zoo)
library(MASS)
library(plotly)
library(pracma)


N = 37
nI = 100

reynolds = 50

pRelax = 0.5
uRelax = 0.5
vRelax = 0.5

d = N / reynolds

P = matrix(data = 0, N, N)
U <- list()
U$face = matrix(data = 0, N , N - 1)
V <- list()
V$face = matrix(data = 0, N - 1, N)

U$bottom = 1
U$left = 0
U$right = 0
U$top = 0

V$bottom = 0
V$left = 0
V$right = 0
V$top = 0

hybridScheme <- function(vn,uw,ue,vs,dn,dw,de,ds){
  aCoeffs <- list()
  aCoeffs$n = max(-vn,dn-vn/2,0)
  aCoeffs$w = max(uw,dw+uw/2,0)
  aCoeffs$e = max(-ue,de-ue/2,0)
  aCoeffs$s = max(vs,ds+vs/2,0)
  return(aCoeffs)
}

upwindScheme <- function(vn,uw,ue,vs,dn,dw,de,ds){
  aCoeffs <- list()
  aCoeffs$n = dn + max(0, -vn)
  aCoeffs$w = dw + max(0, uw)
  aCoeffs$e = de + max(0, -ue)
  aCoeffs$s = ds + max(0, vs)
  return(aCoeffs)
}

uSolve <- function(uFace, vFace, p, uFaceHorizontalBounded, vFaceVerticalBounded) {
  uNode = apply(uFaceHorizontalBounded, 1, function(x)
    rollmean(x, 2))
  vCorner = apply(vFaceVerticalBounded, 1, function(x)
    rollmean(x, 2))
  
  vCornerLinear = c(vCorner)
  uNodeLinear = c(uNode)
  
  A <- matrix(0, N * (N - 1), N * (N - 1) + 2 * N - 2)
  
  for (i in 1:(N * (N - 1))) {
    Dn = d
    Dw = d
    De = d
    Ds = d
    l = (i - 1) %/% (N - 1)
    Vn = vCornerLinear[i]
    Uw = uNodeLinear[i+l]
    Ue = uNodeLinear[i + 1 + l]
    Vs = vCornerLinear[i + N - 1]
    bcn = 1
    bcw = 1
    bce = 1
    bcs = 1
    if (i <= N - 1) {
      Dn = 2 * d
      bcn = U$top
    }
    if ((i - 1) %% (N - 1) == 0) {
      bcw = U$left
    }
    if (i %% (N - 1) == 0) {
      bce = U$right
    }
    if (i > (N - 1) * (N - 1)) {
      Ds = 2 * d
      bcs = U$bottom
    }
    a = hybridScheme(Vn,Uw,Ue,Vs,Dn,Dw,De,Ds)
    
    ap = -(a$n + a$w + a$e + a$s + Ue - Uw + Vn - Vs)/uRelax
    
    A[i, i] =  a$n * bcn
    A[i, i + N - 2] = a$w * bcw
    A[i, i + N - 1] = ap
    A[i, i + N] = a$e * bce
    A[i, i + 2 * N - 2] = a$s * bcs
  }
  squareRange = N:(N + N * (N - 1) - 1)
  AU = A[, squareRange]
  BU = -rowSums(A[, -(squareRange)])
  ApU = t(matrix(-diag(AU), N - 1, N))
  BU = BU + c(apply(p, 1, diff)) - (1-uRelax)*c(t(ApU*uFace))
  uLinear = solve(AU, BU)
  uFaceNew = t(matrix(uLinear, N - 1, N))
  return(list(uFaceNew, ApU))
}

vSolve <- function(uFace, vFace, p, uFaceHorizontalBounded, vFaceVerticalBounded) {
  uCorner = apply(uFaceHorizontalBounded, 2, function(x)
    rollmean(x, 2))
  vNode = apply(vFaceVerticalBounded, 2, function(x)
    rollmean(x, 2))
  
  uCornerLinear = c(t(uCorner))
  vNodeLinear = c(t(vNode))
  
  A <- matrix(0, N * (N - 1), N * (N - 1) + 2 * N)
  
  for (i in 1:(N * (N - 1))) {
    Dn = d
    Dw = d
    De = d
    Ds = d
    l = (i - 1) %/% N
    Vn = vNodeLinear[i]
    Uw = uCornerLinear[i + l]
    Ue = uCornerLinear[i + 1 + l]
    Vs = vNodeLinear[i + N]
    bcn = 1
    bcw = 1
    bce = 1
    bcs = 1
    if (i <= N) {
      bcn = V$top
    }
    if ((i - 1) %% N == 0) {
      Dw = 2 * d
      bcw = V$left
    }
    if (i %% N == 0) {
      De = 2 * d
      bce = V$right
    }
    if (i > N*(N-2)) {
      bcs = V$bottom
    }
    a = hybridScheme(Vn,Uw,Ue,Vs,Dn,Dw,De,Ds)
    ap = -(a$n + a$w + a$e + a$s + Ue - Uw + Vn - Vs)/vRelax
    
    A[i, i] =  a$n * bcn
    A[i, i + N - 1] = a$w * bcw
    A[i, i + N ] = ap
    A[i, i + N + 1] = a$e * bce
    A[i, i + 2 * N] = a$s * bcs
  }
  squareRange = (N+1):(N + N * (N - 1))
  AV = A[, squareRange]
  ApV = t(matrix(-diag(AV), N, N - 1))
  BV = -rowSums(A[, -(squareRange)])
  BV =  BV  + c(t(apply(p, 2, diff))) - (1-vRelax)*c(t(ApV*vFace))
  vLinear = solve(AV, BV)
  vFaceNew = t(matrix(vLinear, N, N - 1))
  return(list(vFaceNew, ApV))
}

pSolve <- function(apU, apV, uFace, vFace) {
  apU = apU[2:(N - 1),]
  apU = c(t(apU))
  apV = apV[, 2:(N - 1)]
  apV = c(t(apV))
  uFaceStripped = uFace[2:(N - 1),]
  vFaceStripped = vFace[,2:(N-1)]
  a = 1 / apU
  b = 1 / apV
  A <- matrix(0, (N - 2) * (N - 2), (N - 2) * (N - 2) + 2*(N-2))
  for (i in 1:((N - 2) * (N - 2))) {
    l = ((i - 1) %/% (N - 2))
    an = b[i]
    aw = a[i + l]
    ae = a[i + 1 + l]
    as = b[i + N - 2]
    ap = -(an + aw + ae + as)
    if (i <= N - 2) {
      ap =  ap + an
      an = 0
    }
    if ((i - 1) %% (N - 2) == 0) {
      ap = ap + aw
      aw = 0
    }
    if (i %% (N - 2) == 0) {
      ap = ap + ae
      ae = 0
    }
    if (i > (N - 2) * (N - 3)) {
      ap = ap + as
      as = 0
    }
    A[i, i] = an
    A[i, i + N - 3] = aw
    A[i, i + N - 2] = ap
    A[i, i + N  -1] = ae
    A[i, i + 2 * N - 4] = as
    
  }
  squareRange = (N-1):(N-2+(N-2)*(N-2))
  A = A[, squareRange]
  A = rbind(A,1)
  dU = c(apply(uFaceStripped, 1, diff))
  dV = c(t(apply(vFaceStripped, 2, diff)))
  B = dU + dV
  B = c(B,0)
  P = ginv(A) %*% B
  P = t(matrix(P, N - 2, N - 2))
  P = rbind(P[1, ], P, P[N - 2, ])
  P = cbind(P[, 1], P, P[, N - 2])
  P[1, 1] = (P[1, 2] + P[2, 1]) / 2
  P[N, 1] = (P[N - 1, 1] + P[N, 2]) / 2
  P[1, N] = (P[2, N] + P[1, N - 1]) / 2
  P[N, N] = (P[N - 1, N] + P[N, N - 1]) / 2
  return(P)
}

plotVelocityField <- function(){
  U$plot = cbind(U$left,U$face,U$right)
  U$plot = rbind(U$top,U$plot,U$bottom)
  U$plot = t(apply(U$plot, 1, function(x)
    rollmean(x, 2)))
  U$plot = cbind(U$left,U$plot,U$right)
  
  V$plot = rbind(V$top,V$face,V$bottom)
  V$plot = cbind(V$left,V$plot,V$right)
  V$plot = apply(V$plot, 2, function(x)
    rollmean(x, 2))
  V$plot = rbind(V$top,V$plot,V$bottom)
  plot(range(X),range(X),type="n")
  quiver(G$X,G$Y,U$plot,V$plot,scale=0.5,length=0.1,angle=10)
}

Xp = seq(1/(2*N),1-1/(2*N),1/N)
X = c(0,Xp,1)
G <- meshgrid(X, X)

x11()

for (i in 1:nI) {
  U$FaceHorizontalBounded = cbind(U$left, U$face, U$right)
  V$FaceVerticalBounded = rbind(V$top, V$face, V$bottom)
  U$solution = uSolve(U$face, V$face, P, U$FaceHorizontalBounded, V$FaceVerticalBounded)
  U$faceNew = U$solution[[1]]
  U$Ap = U$solution[[2]]
  V$solution = vSolve(U$face, V$face, P, U$FaceHorizontalBounded, V$FaceVerticalBounded)
  V$faceNew = V$solution[[1]]
  V$Ap = V$solution[[2]]
  pline = pSolve(U$Ap, V$Ap, U$faceNew, V$faceNew)
  P = P + pRelax*pline
  U$face = U$faceNew - t(apply(pline, 1, diff))/U$Ap
  V$face = V$faceNew - apply(pline, 2, diff)/V$Ap
  plotVelocityField()
}

plot_ly(x = Xp, y = Xp, z = P) %>% add_surface()

#plot_ly(x = X, y = X, z = U$plot) %>% add_surface()
#plot_ly(x = X, y = X, z = V$plot) %>% add_surface()
