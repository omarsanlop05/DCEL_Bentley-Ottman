# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 2026

@author: AlexL
"""

from Punto import Punto

class Linea:
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def distancia(self, punto):
        Ap = -self.B
        Bp = self.A
        Cp = self.B*punto.x - self.A*punto.y

        Lp = Linea(Ap,Bp,Cp)

        inter = self.interseccion(Lp)

        return inter.distancia(punto)
    
    def interseccion(self,linea):
        n = (self.A * linea.B) - (linea.A * self.B)
        if abs(n) < 1e-9:
            return None

        detX = -(self.C * linea.B) + (linea.C * self.B)
        detY = -(self.A * linea.C) + (linea.A * self.C)

        x = detX / n
        y = detY / n

        return Punto(x,y)