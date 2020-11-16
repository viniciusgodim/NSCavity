import numpy as np
from numpy.linalg import solve,lstsq
from scipy.ndimage.filters import uniform_filter1d
import matplotlib.pyplot as plt

N = 50

vRange = np.array((3+N)*[False]+(N*[True]+2*[False])*(N-1)+(N+1)*[False])
uRange = np.array((3+N-1)*[False]+((N-1)*[True]+2*[False])*N+N*[False])
pRange = np.array((3+N-2)*[False]+((N-2)*[True]+2*[False])*(N-2)+(N-1)*[False])

Xp = np.linspace(1/(2*N),1-1/(2*N),N)
X = np.append(np.append(0,Xp),1)

nIterations = 100

reynolds = 25

pRelax = 0.5
uRelax = 0.5
vRelax = 0.5

d = N/reynolds

P = np.zeros((N,N))
Uface = np.zeros((N,N-1))
Vface = np.zeros((N-1,N))

uBC = [0,0,0,1]
uBCArray = [np.ones((1,N+2))*uBC[0],np.ones((N,1))*uBC[1],np.ones((N,1))*uBC[2],np.ones((1,N+2))*uBC[3]]

vBC = [0,1,0,0]
vBCArray = [np.ones((1,N))*vBC[0],np.ones((N+2,1))*vBC[1],np.ones((N+2,1))*vBC[2],np.ones((1,N))*vBC[3]]

def movingAverage(x):
    return np.convolve(x, np.ones(2)/2, mode='valid')

def hybridScheme(f,d):
    return np.array([max(f[i],d[i]+f[i]/2,0) for i in range(4)])

def uSolve(uFace,vFace,p):
    uNode = np.apply_along_axis(movingAverage, 1, uFaceHorizontalBounded)
    vCorner =  np.apply_along_axis(movingAverage, 1, vFaceVerticalBounded)
    uNodeLinear = uNode.flatten()
    vCornerLinear = vCorner.flatten()
    A = np.zeros((N * (N - 1), (N+2)*(N+1)))
    for i in range(N*(N-1)):
        D = np.ones(4)* d
        l = i // (N-1)
        F = [-vCornerLinear[i],uNodeLinear[i+l],-uNodeLinear[i+1+l],vCornerLinear[i + N - 1]]
        wBC = np.ones(4)
        if i + 1 <= N - 1:
            D[0] *= 2
            wBC[0] = uBC[0]
        if i % (N-1) == 0:
            wBC[1] = uBC[1]
        if (i+1) % (N-1) == 0:
            wBC[2] = uBC[2]
        if i + 1 > (N-1)*(N-1):
            D[3] *= 2
            wBC[3] = uBC[3]
        a = hybridScheme(F,D)
        ap = -(sum(a) - sum(F))/uRelax
        aW = wBC * a
        A[i, i + 1 + 2*l] =  aW[0]
        A[i, i + 1 + 2*l + N] = aW[1]
        A[i, i + 2 + 2*l + N] = ap
        A[i, i + 3 + 2*l + N] = aW[2]
        A[i, i + 3 + 2*l + 2*N] = aW[3]
    AU = A[:,uRange]
    BU = A[:,~uRange]
    BU = -np.sum(BU,axis=1)
    ApU = -np.diagonal(AU)
    ApU = np.reshape(ApU,(-1,N-1))
    BU = BU + np.diff(p,axis=1).flatten() - (1-uRelax)*(ApU*uFace).flatten()
    uLinear = solve(AU,BU)
    uFaceNew = np.reshape(uLinear,(-1,N-1))
    return [uFaceNew,ApU]

def vSolve(uFace,vFace,p):
    uCorner = np.apply_along_axis(movingAverage, 0, uFaceHorizontalBounded)
    vNode =  np.apply_along_axis(movingAverage, 0, vFaceVerticalBounded)
    uCornerLinear = uCorner.flatten()
    vNodeLinear = vNode.flatten()
    A = np.zeros((N * (N - 1), (N+2) * (N+1)))
    for i in range(N*(N-1)):
        D = np.ones(4)* d
        l = i // N
        F = [-vNodeLinear[i],uCornerLinear[i+l],-uCornerLinear[i+1+l],vNodeLinear[i + N]]
        wBC = np.ones(4)
        if i + 1 <= N:
            wBC[0] = vBC[0]
        if i % N == 0:
            D[1] *= 2
            wBC[1] = vBC[1]
        if (i+1) % N == 0:
            D[2] *= 2
            wBC[2] = vBC[2]
        if (i + 1) > (N-2)*N:
            wBC[3] = vBC[3]
        a = hybridScheme(F,D)
        ap = -(sum(a) - sum(F))/vRelax
        aW = wBC * a
        A[i, i + 1 + 2*l] =  aW[0]
        A[i, i + 1 + 2*l + N + 1] = aW[1]
        A[i, i + 2 + 2*l + N + 1] = ap
        A[i, i + 3 + 2*l + N + 1] = aW[2]
        A[i, i + 5 + 2*l + 2*N] = aW[3]
    AV = A[:,vRange]
    BV = A[:,~vRange]
    BV = -np.sum(BV,axis=1)
    ApV = -np.diagonal(AV)
    ApV = np.reshape(ApV,(-1,N))
    BV = BV + np.diff(p,axis=0).flatten() - (1-vRelax)*(ApV*vFace).flatten()
    vLinear = solve(AV,BV)
    vFaceNew = np.reshape(vLinear,(-1,N))
    return [vFaceNew,ApV]

def pSolve(apU,apV,uFace,vFace):
    apU = apU[range(1,N-1),:].flatten()
    apV = apV[:,range(1,N-1)].flatten()
    uFaceStripped = uFace[range(1,N-1),:]
    vFaceStripped = vFace[:,range(1,N-1)]
    q = 1/apU
    r = 1/apV
    A = np.zeros(((N - 2)**2, N**2))
    for i in range((N-2)**2):
        l = (i // (N - 2))
        a = [r[i],q[i+l],q[i+1+l],r[i+N-2]]
        ap = -sum(a)
        conditions = [i + 1 <= N - 2,i % (N-2) == 0,(i+1)%(N-2)==0,(i+1)>(N-2)*(N-3)]
        for k,cond in enumerate(conditions):
            if cond:
                ap = ap + a[k]
                a[k] = 0
        A[i, i + 1 + 2*l] =  a[0]
        A[i, i + 1 + 2*l + N - 1] = a[1]
        A[i, i + 2 + 2*l + N - 1] = ap
        A[i, i + 3 + 2*l + N - 1] = a[2]
        A[i, i + 1 + 2*l + 2*N] = a[3]
    A = A[:,pRange]
    A = np.vstack((A,np.ones((1,(N-2)*(N-2)))))
    dU = np.diff(uFaceStripped,axis=1).flatten()
    dV = np.diff(vFaceStripped,axis=0).flatten()
    B = dU + dV
    B = np.append(B,0)
    P = lstsq(A,B,rcond=None)[0]
    P = np.reshape(P,(-1,N-2))
    P = np.vstack(([P[0,:]],P,[P[-1,:]]))
    P = np.hstack(([[pr] for pr in P[:,0]],P,[[pr] for pr in P[:,-1]]))
    P[0, 0] = (P[0, 1] + P[1, 0]) / 2
    P[-1, 0] = (P[-2, 0] + P[-1, 1]) / 2
    P[0,-1] = (P[1,-1] + P[0,-2]) / 2
    P[-1, -1] = (P[-2,-1] + P[-1,-2]) / 2
    return P

for i in range(nIterations):
    print(i)
    uFaceHorizontalBounded = np.hstack((uBCArray[1],Uface,uBCArray[2]))
    vFaceVerticalBounded = np.vstack((vBCArray[0],Vface,vBCArray[3]))
    USolution = uSolve(Uface,Vface,P)
    UfaceNew = USolution[0]
    UAp = USolution[1]
    VSolution = vSolve(Uface,Vface,P)
    VfaceNew = VSolution[0]
    VAp = VSolution[1]
    pline = pSolve(UAp,VAp,UfaceNew,VfaceNew)
    P = P + pRelax*pline
    Uface = UfaceNew - np.diff(pline,axis=1)/UAp
    Vface = VfaceNew - np.diff(pline,axis=0)/VAp

Uplot = np.hstack((uBCArray[1],Uface,uBCArray[2]))
Uplot = np.apply_along_axis(movingAverage, 1, Uplot)
Uplot = np.hstack((uBCArray[1],Uplot,uBCArray[2]))
Uplot = np.vstack((uBCArray[0],Uplot,uBCArray[3]))
Vplot = np.vstack((vBCArray[0],Vface,vBCArray[3]))
Vplot = np.apply_along_axis(movingAverage, 0, Vplot)
Vplot = np.vstack((vBCArray[0],Vplot,vBCArray[3]))
Vplot = np.hstack((vBCArray[1],Vplot,vBCArray[2]))
fig, ax = plt.subplots()
ax.quiver(X, X, Uplot, Vplot)
plt.show()

XP,YP = np.meshgrid(Xp,Xp)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(XP, YP, P)
plt.show()
