import numpy as np
from matplotlib import pyplot as plt

def f(x):
    res = 0
    for i in range(0, 19):
        if res > 9:
            res = res - i*x[i]*x[i];
        else:
            res = res + i*x[i]*x[i];
    return res

itor = 200
dim = 20
xlist = np.zeros((itor, dim))
x = np.random.randn(dim)
step = 100
for i in range(itor):
    dir = np.random.randn(dim)
    cnt = 0
    dir = dir / np.linalg.norm(dir, ord=1, axis=0)
    step = 0.95 * step
    d0 = step / 2
    #print(step)
    xl = x + d0 * dir / 2
    xr = x - d0 * dir / 2
    fleft = f(xl)
    fright = f(xr)
    tx = x - step * dir * np.sign(fleft - fright)
    xlist[i, :] = tx
    if abs(tx[0]) > 100 or abs(tx[1]) > 100:
        continue
    else:
        x = tx
    #print(xl)
    #print(min(fleft, fright))
    #print(x)
    #print(f(x))

if fleft < fright:
    print(fleft)
else:
    print(fright)
print(x)
print(f(x))
plt.figure(0)
plt.plot(xlist[:,0], xlist[:,1])
plt.show()
