# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 2026

@author: AlexL
"""

import numpy as np
from Linea import Linea

class Segmento:
    def __init__(self, p1, p2):
        if p1.y < p2.y:
            self.p1, self.p2 = p2, p1
        elif p1.y == p2.y and p1.x > p2.x:
            self.p1, self.p2 = p2, p1
        else:
            self.p1, self.p2 = p1, p2

    def __hash__(self):
        return hash((self.p1, self.p2))

    def __eq__(self, other):
        if not isinstance(other, Segmento):
            return NotImplemented
        return (self.p1 == other.p1 and self.p2 == other.p2) or \
            (self.p1 == other.p2 and self.p2 == other.p1)
    
    def longitud(self):
        x = (self.p1.x - self.p2.x)**2
        y = (self.p1.y - self.p2.y)**2
        return np.sqrt(x+y)

    def pendiente(self):
        return (self.p2.y - self.p1.y) / (self.p2.x - self.p1.x)

    def aLinea(self):
        A = (self.p1.y - self.p2.y)
        B = (self.p2.x - self.p1.x)
        C = (self.p1.x*self.p2.y) - (self.p1.y*self.p2.x)

        return Linea(A,B,C)
    
    def puntoEnSegmento(self, p):
        return (min(self.p1.x, self.p2.x) <= p.x <= max(self.p1.x, self.p2.x) and
                min(self.p1.y, self.p2.y) <= p.y <= max(self.p1.y, self.p2.y))

    def colinealCon(self, seg):
        # área del triángulo = 0 → colineales
        return ((seg.p1.x - self.p1.x) * (self.p2.y - self.p1.y) -
                (seg.p1.y - self.p1.y) * (self.p2.x - self.p1.x)) == 0
    
    def distancia(self, punto):
        l = self.aLinea()

        Ap = -l.B
        Bp = l.A
        Cp = l.B*punto.x - l.A*punto.y

        Lp = Linea(Ap,Bp,Cp)

        inter = l.interseccion(Lp)

        dentro = self.puntoEnSegmento(inter)

        if dentro:
            return inter.distancia(punto), inter
        else:
            # Regresar el punto (p1 o p2) más cercano al punto

            d1 = self.p1.distancia(punto)
            d2 = self.p2.distancia(punto)

            if d1 < d2:
                return d1, self.p1
            else:
                return d2, self.p2
    
    ## Convertir a linea, después revisar si la intersección está entre los segmentos
    def interseccion(self,segmento):
        l = self.aLinea()
        ls = segmento.aLinea()

        p = l.interseccion(ls)

        if p is not None:
            if self.puntoEnSegmento(p) and segmento.puntoEnSegmento(p):
                return [p]
            else:
                return []

        if not self.colinealCon(segmento):
            return []

        pts = []
        for pto in [self.p1, self.p2, segmento.p1, segmento.p2]:
            if self.puntoEnSegmento(pto) and segmento.puntoEnSegmento(pto):
                pts.append(pto)

        return pts





        