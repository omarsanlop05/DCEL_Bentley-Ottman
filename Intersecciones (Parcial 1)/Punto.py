# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 14:37:28 2026

@author: AlexL
"""

import numpy as np

class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        if not isinstance(other, Punto):
            return NotImplemented
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __lt__(self, other):
        if self.y != other.y:
            return self.y > other.y  # Y descendente (de arriba hacia abajo)
        return self.x < other.x  # desempate por X ascendente

    def rotar(self, alpha):
        x_nuevo = self.x * np.cos(alpha) - self.y * np.sin(alpha)
        y_nuevo = self.x * np.sin(alpha) + self.y * np.cos(alpha)
        self.x = x_nuevo
        self.y = y_nuevo

    def trasladar(self, dx, dy):
        self.x += dx
        self.y += dy

    def comparar(self, punto):
        return self.x == punto.x and self.y == punto.y
    
    def aPolar(self):
        r = np.sqrt(self.x**2 + self.y**2)
        theta = np.arctan2(self.y, self.x)

        return r, theta
    
    def distancia(self, otro):
        x = (self.x - otro.x)**2
        y = (self.y - otro.y)**2

        d = np.sqrt(x+y)

        return d

    def __str__(self):
        return f"({self.x:.3f}, {self.y:.3f})"
    

