import numpy as np

def findSpan(n,p,u,U):
    if (u == U[n+1]):
        return n
    low = p
    high = n+1
    mid = (low+high)/2.0
    while (u<U[mid] || u>U[mid+1]):
        if (u<U[mid]):
            high = mid
        else:
            low = mid
        mid = np.floor((low+high)/2.0)
    return mid
