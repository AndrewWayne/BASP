import numpy as np
def  zuo(a,b):
    return (1+(a+b+1)**2*(19-14*a+3*b**2-14*b+6*a*b+3*b**2))*(30+(2*a-3*b)**2*(18-32*a+12*a**2+48*b-36*a*b+27*b**2))
def  you(c,d):
    return(1+(c+d+1)**2*(19-14*c+3*d**2-14*d+6*c*d+3*d**2))*(30+(2*c-3*d)**2*(18-32*c+12*c**2+48*d-36*c*d+27*d**2))
x=np.random.randn(2)
step=2
for i in range (500):
    dir = np.random.randn(2)
    dir = dir / np.linalg.norm(dir, ord=1, axis=0)
    step =0.95*step
    d0=step/2
    xl=x+d0*dir/2
    xr=x-d0*dir/2
    a=xl[0]
    b=xl[1]
    c=xr[0]
    d=xr[1]
    fleft=zuo(a,b)
    fright=you(c,d)
    x=x-step*dir*np.sign(fleft-fright)
    if abs(x[0])>2 or abs(x[1])>2:
        continue
if fleft<fright:
    print(fleft)
else:
    print(fright)
print(x)
