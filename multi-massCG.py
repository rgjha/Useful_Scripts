# Written by - Vasudev Shyam 
# and slightly edited by me.
# Jan 9, 2020
# Reference: Page 6-7 of https://arxiv.org/pdf/hep-lat/9612014.pdf

import numpy as np
from numpy import linalg as LA

def dagger(a):
    return np.transpose(a).conj()

def mmconjgrad(A, b, s):
    n = len(b)
    l = len(s)
    xs = np.zeros((l,n), dtype = float)
    ps = np.zeros((l,n), dtype = float)
    alpha_s = np.ones(l, dtype = float)
    beta_s = np.ones(l, dtype = float)
    zeta_k = np.ones(l, dtype = float)
    zeta_km = np.ones(l, dtype = float)
    zeta_kp = np.ones(l, dtype = float)
    co = np.ones(l, dtype = float)
    ct = np.ones(l, dtype = float)


    zet_k = zet_km = zet_kp = 1
    bet = 0
    alpha_km = 1

    for j in range(l):
        ps[j] = b
        beta_s[j] = bet
        alpha_s[j] = alpha_km
        zeta_k[j] = zet_k
        zeta_km[j] = zet_km
        zeta_kp[j] = zet_kp

    r = b
    p = b

    for i in range(l):
        for j in range(100):
            ap = np.dot(A,p)
            r_k_norm = np.dot(r, r)
            alpha_k = -r_k_norm / np.dot(p, ap)
            zeta_kp[i] = zeta_k[i] * zeta_km[i] * alpha_km
            co[i] = alpha_k * bet * (zeta_km[i] - zeta_k[i])
            ct[i] = zeta_km[i] * alpha_km * (1 - s[i] * alpha_k)
            zeta_kp[i] /= co[i] + ct[i]
            alpha_s[i] = (alpha_k * zeta_kp[i]) / zeta_k[i]

            xs[i] -= alpha_s[i] * ps[i]
            # xs[i] is psim[j] in SUSY code 
            # https://github.com/daschaich/susy/blob/master/4d_Q16/susy/congrad_multi.c
            r += alpha_k * ap
            r_kp_norm = np.dot(r, r)

            bet = r_kp_norm / r_k_norm

            beta_s[i] = bet*(zeta_kp[i]/zeta_k[i])*(alpha_s[i]/alpha_k)

            if r_kp_norm < 1e-15:
                print("Finishing up at Iteration", j)
                zet_k = zet_km = zet_kp = 1
                r = ps[i-l+1]
                p = ps[i-l+1]
                bet = 0
                alpha_km = 1
                break
            else:
                # Scroll the variables for next iteration
                alpha_km = alpha_k
                zeta_km[i] = zeta_k[i]
                zeta_k[i] = zeta_kp[i]
            ps[i] = (r*zeta_kp[i])+(beta_s[i]*p)

    return xs


R = np.random.rand(5,5)
A = (R + R.T) * 0.5
k = 10.0
A += k * np.identity(5)
b = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
s = np.array([0.5*k,0.6*k, 0.4*k, 0.7*k])
B0 = A + np.dot(s[0],np.eye(5))
B1 = A + np.dot(s[1],np.eye(5))
B2 = A + np.dot(s[2],np.eye(5))
B3 = A + np.dot(s[3],np.eye(5))

print(np.dot(LA.inv(B0),b))
print(np.dot(LA.inv(B1),b))
print(np.dot(LA.inv(B2),b))
print(np.dot(LA.inv(B3),b))

m = mmconjgrad(A,b,s)
print(m)

#print ("Is CG result = inversion result-> ", np.allclose(np.dot(np.linalg.inv(B),b), result))


