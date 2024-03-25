import numpy as np
import numba as nb

def make_squares(rsqr, coeffs):
    n = len(coeffs)
    res = 0.0
    for i in range(n):
        res += rsqr**(i+1)*coeffs[i]
    return res

class DistorsionPolynome(object):
    def __init__(self, k_coefficients,p_coefficients,x_c=0,y_c=0):
        self.k_coefficients = k_coefficients
        self.p_coefficients = p_coefficients
        self.x_c = x_c
        self.y_c = y_c

    def __tangent_rterm(self,rsqr):
        if len(self.p_coefficients)==2:
            return 1.0
        elif len(self.p_coefficients)>2:
            return 1.0 + make_squares(rsqr,self.p_coefficients[2:])
        else:
            return 0.0

    def __p1_p2(self):
        if self.p_coefficients:
            p1 = self.p_coefficients[0]
        else:
            p1 = 0.0

        if len(self.p_coefficients)>1:
            p2 = self.p_coefficients[1]
        else:
            p2 = 0.0

        return p1,p2

    def undistort(self,x_d,y_d):
        rsqr = x_d*x_d + y_d*y_d
        dx = x_d - self.x_c
        dy = y_d - self.y_c
        rterm = self.__tangent_rterm(rsqr)
        p1,p2 = self.__p1_p2()
        x_u = x_d + dx*(make_squares(rsqr,self.k_coefficients)) + (p1*(rsqr + 2*dx**2) + 2*p2*dx*dy)*rterm
        y_u = y_d + dy*(make_squares(rsqr,self.k_coefficients)) + (2*p1*dx*dy + p2*(rsqr + 2*dy**2))*rterm
        return x_u, y_u
