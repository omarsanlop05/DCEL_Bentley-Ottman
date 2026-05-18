import numpy as np
import math

class Arista:
    def __init__(self, nombre):
        self.nombre = nombre
        self.origen = None
        self.antiarista = None
        self.cara = None
        self.siguiente = None
        self.anterior = None

    def __str__(self):

        origen = self.origen.nombre if self.origen else "None"
        pareja = self.antiarista.nombre if self.antiarista else "None"
        cara = self.cara.nombre if self.cara else "None"
        sig = self.siguiente.nombre if self.siguiente else "None"
        ant = self.anterior.nombre if self.anterior else "None"

        return (f"Arista {self.nombre}: "
                f"origen={origen}, pareja={pareja}, "
                f"cara={cara}, siguiente={sig}, anterior={ant}")

class Cara:
    def __init__(self, nombre):
        self.nombre = nombre
        self.aristas_internas = []
        self.aristas_externas = None

    def __str__(self):
        aristas_externas = (
            self.aristas_externas.nombre
            if self.aristas_externas else
            "None"
        )
        aristas_internas = [
            e.aristas_internas.nombre
            for e in self.aristas_internas
        ]
        return (
            f"{self.nombre} "
            f"{aristas_externas} "
            f"{aristas_internas}"
        )

class Vertice:
    def __init__(self, nombre, pt):
        self.nombre = nombre
        self.pt = pt
        self.arista_adyacente = None

    def __str__(self):
        inc = self.arista_adyacente.nombre if self.arista_adyacente else "None"
        return f"Vertice {self.nombre}: punto={self.pt}, incidente={inc}"


class Punto:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, otro):
        if not isinstance(otro, Punto): return False
        return math.isclose(self.x, otro.x, abs_tol=1e-7) and \
            math.isclose(self.y, otro.y, abs_tol=1e-7)

    def __hash__(self):
        # Redondeamos para que puntos "casi iguales" generen el mismo hash
        return hash((round(self.x, 7), round(self.y, 7)))

    def __lt__(self, other):
        if self.y != other.y:
            return self.y > other.y  # Y descendente (de arriba hacia abajo)
        return self.x < other.x  # desempate por X ascendente

    def comparar(self, otro):
        return self.__eq__(otro)

    def rotar(self, alpha):
        x_nuevo = self.x * np.cos(alpha) - self.y * np.sin(alpha)
        y_nuevo = self.x * np.sin(alpha) + self.y * np.cos(alpha)
        self.x = x_nuevo
        self.y = y_nuevo

    def trasladar(self, dx, dy):
        self.x += dx
        self.y += dy

    def aPolar(self):
        r = np.sqrt(self.x ** 2 + self.y ** 2)
        theta = np.arctan2(self.y, self.x)

        return r, theta

    def distancia(self, otro):
        x = (self.x - otro.x) ** 2
        y = (self.y - otro.y) ** 2

        d = np.sqrt(x + y)

        return d

    def __str__(self):
        return f"({self.x:.3f}, {self.y:.3f})"


class Linea:
    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def distancia(self, punto):
        Ap = -self.B
        Bp = self.A
        Cp = self.B * punto.x - self.A * punto.y

        Lp = Linea(Ap, Bp, Cp)

        inter = self.interseccion(Lp)

        return inter.distancia(punto)

    def interseccion(self, linea):
        n = (self.A * linea.B) - (linea.A * self.B)
        if abs(n) < 1e-12:
            return None

        detX = -(self.C * linea.B) + (linea.C * self.B)
        detY = -(self.A * linea.C) + (linea.A * self.C)

        x = detX / n
        y = detY / n

        return Punto(x, y)


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
        x = (self.p1.x - self.p2.x) ** 2
        y = (self.p1.y - self.p2.y) ** 2
        return np.sqrt(x + y)

    def pendiente(self):
        return (self.p2.y - self.p1.y) / (self.p2.x - self.p1.x)

    def aLinea(self):
        A = (self.p1.y - self.p2.y)
        B = (self.p2.x - self.p1.x)
        C = (self.p1.x * self.p2.y) - (self.p1.y * self.p2.x)

        return Linea(A, B, C)

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
        Cp = l.B * punto.x - l.A * punto.y

        Lp = Linea(Ap, Bp, Cp)

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
    def interseccion(self, segmento):
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

class Evento:
    def __init__(self, punto, inicio = None, contiene = None):
        self.punto = punto

        if inicio is None:
            self.inicio = []
        else:
            self.inicio = inicio

        if contiene is None:
            self.contiene = []
        else:
            self.contiene = contiene

    def __lt__(self, other):
        if self.punto.y != other.punto.y:
            return self.punto.y > other.punto.y
        return self.punto.x < other.punto.x

    def __hash__(self):
        return hash((self.punto.x, self.punto.y))

    def __eq__(self, other):
        if not isinstance(other, Evento):
            return NotImplemented
        return self.punto.x == other.punto.x and self.punto.y == other.punto.y

    def agregarInicial(self, segmento):
        self.inicio.append(segmento)

    def agregarContiene(self, segmento):
        self.contiene.append(segmento)

def leerVertices(ruta, vertices):
    datos_vertices = []

    with open(ruta + ".vertices", "r") as f:
        for linea in f:
            if linea.startswith("#") or linea.strip() == "" or linea.startswith("Nombre") or linea.startswith("Archivo") or linea.startswith("-"):
                continue

            datos = linea.split()
            nombre = datos[0]
            x = float(datos[1])
            y = float(datos[2])
            adyacente = datos[3]

            punto = Punto(x, y)
            v = Vertice(nombre, punto)

            vertices[nombre] = v
            datos_vertices.append((nombre, adyacente.strip()))

    return datos_vertices

def leerAristas(ruta, aristas):
    datos_aristas = []

    with open(ruta + ".aristas", "r") as f:
        for linea in f:
            if linea.startswith("#") or linea.strip() == "" or linea.startswith("Nombre") or linea.startswith("Archivo") or linea.startswith("-"):
                continue

            datos = linea.split()

            nombre = datos[0]
            origen = datos[1]
            antiarista = datos[2]
            cara = datos[3]
            siguiente = datos[4]
            anterior = datos[5]

            a = Arista(nombre)
            aristas[nombre] = a

            datos_aristas.append((
                nombre.strip(),
                origen.strip(),
                antiarista.strip(),
                cara.strip(),
                siguiente.strip(),
                anterior.strip()
            ))

    return datos_aristas

def leerCaras(ruta, caras):
    datos_caras = []

    with open(ruta + ".caras", "r") as f:
        for linea in f:
            if linea.startswith("#") or linea.strip() == "" or linea.startswith("Nombre") or linea.startswith("Archivo") or linea.startswith("-"):
                continue

            datos = linea.split()
            nombre = datos[0]
            internas = None if datos[1] == "None" else datos[1]
            externas = None if datos[2] == "None" else datos[2]

            c = Cara(nombre.strip())
            caras[nombre.strip()] = c

            datos_caras.append((nombre.strip(), internas, externas))

    return datos_caras

def definirObjetos(ruta):

    vertices = {}
    aristas = {}
    caras = {}

    data_aristas = leerAristas(ruta, aristas)
    data_vertices = leerVertices(ruta, vertices)
    data_caras = leerCaras(ruta, caras)

    # Aristas
    for nombre, origen, antiarista, cara, siguiente, anterior in data_aristas:

        a = aristas[nombre]
        a.origen = vertices.get(origen)
        a.antiarista = aristas.get(antiarista)
        a.cara = caras.get(cara)
        a.siguiente = aristas.get(siguiente)
        a.anterior = aristas.get(anterior)

    # Aristas
    for nombre, antiarista in data_vertices:
        if antiarista in aristas:
            vertices[nombre].arista_adyacente = aristas[antiarista]
        else:
            print(f"⚠️ Arista incidente '{antiarista}' no existe en vértice {nombre}")

    # Caras
    for nombre, aristas_internas, aristas_externas in data_caras:

        c = caras[nombre]
        if aristas_internas != None:
            aristas_internas = aristas_internas.replace("[", "").replace("]", "").strip()
            lista = [x.strip() for x in aristas_internas.split(",") if x.strip()]

            for ar in lista:
                if ar in aristas:
                    c.aristas_internas.append(aristas[ar])
                else:
                    print(f"⚠️ Arista '{ar}' no existe en cara {nombre}")

        if aristas_externas != None:
            if aristas_externas in aristas:
                c.aristas_externas = aristas[aristas_externas]
            else:
                print(f"⚠️ Arista externa '{aristas_externas}' no existe en cara {nombre}")

    return vertices, aristas, caras