from Figuras import Cara, Punto
import math

class ConstructorCaras:

    def __init__(self, vertices_g, aristas_g):
        self.vertices_g = vertices_g
        self.aristas_g  = aristas_g
        self.caras_g    = {}

    def construir(self):
        ciclos =  self._extraer_ciclos()
        exteriores, agujeros =  self._clasificarCiclos(ciclos)
        contenedor =  self._construirGrafo(agujeros, ciclos)
        self.crearCaras(exteriores,agujeros,contenedor)
        return self.caras_g

    def _extraer_ciclos(self):
        visitadas = set()
        ciclos = []

        for aristas in self.aristas_g.values():
            if aristas.nombre in visitadas:
                continue
            ciclo = []
            actual = aristas

            while actual.nombre not in visitadas:
                ciclo.append(actual)
                visitadas.add(actual.nombre)
                actual = actual.siguiente

            if ciclo:
                ciclos.append(ciclo)

        return ciclos

    def _areaConSigno(self, ciclo):
        area = 0
        for arista in ciclo:
            p1 = arista.origen.pt
            p2 = arista.antiarista.origen.pt
            area += (p1.x * p2.y) - (p2.x * p1.y)
        return area/2

    def _clasificarCiclos(self, ciclos):
        exteriores = []
        agujeros = []

        for ciclo in ciclos:
            if self._areaConSigno(ciclo) > 0:
                exteriores.append(ciclo)
            else:
                agujeros.append(ciclo)
        return exteriores, agujeros


    def verticeMasIzq(self, ciclo):
        mejor = ciclo[0].origen

        for arista in ciclo:
            v = arista.origen
            if v.pt.x < mejor.pt.x:
                mejor = v
            elif v.pt.x == mejor.pt.x and v.pt.y < mejor.pt.y:
                mejor = v
        return mejor

    def _halfEdgeIzq(self, v):
        mejorArista = None
        mejorX = float('-inf')
        y_objetivo = v.pt.y

        for arista in self.aristas_g.values():
            p1 = arista.origen.pt
            p2 = arista.antiarista.origen.pt

            y_min = min(p1.y, p2.y)
            y_max = max(p1.y, p2.y)

            if y_min > y_objetivo or y_max < y_objetivo:
                continue

            if abs(p2.y - p1.y) < 1e-9:
                continue

            t = (y_objetivo - p1.y) / (p2.y -p1.y)
            xCruce =  p1.x + t * (p2.x - p1.x)

            if xCruce >= v.pt.x:
                continue
            if xCruce > mejorX:
                mejorX = xCruce
                mejorArista = arista
        return mejorArista

    def _construirGrafo(self, agujeros, ciclos):
        aristaACiclo = {}
        for ciclo in ciclos:
            for arista in ciclo:
                aristaACiclo[arista.nombre] = ciclo

        contenedor = {}

        for agujero in agujeros:
            vIzq = self.verticeMasIzq(agujero)
            mAIzq = self._halfEdgeIzq(vIzq)

            if mAIzq is None:
                contenedor[id(agujero)] = None
            else:
                contenedor[id(agujero)] = aristaACiclo.get(mAIzq.nombre)
        return contenedor
    def crearCaras(self, exteriores, agujeros, contenedor):
        ##cara infinita
        caraInf = Cara("f_infinita")
        caraInf.aristas_externas = None
        self.caras_g["f_infinita"] = caraInf

        cicloACara = {}

        for i, ciclo in enumerate(exteriores):
            f = Cara(f"f{i}")
            f.aristas_externas = ciclo[0]
            f.aristas_internas = []
            self.caras_g[f.nombre] = f

            for arista in ciclo:
                arista.cara = f

            cicloACara[id(ciclo)] = f

        for agujero in agujeros:
            cicloContenedor = contenedor[id(agujero)]

            if cicloContenedor is None:
                caraContenedor = caraInf

            else:
                caraContenedor = cicloACara.get(id(cicloContenedor), caraInf)

            caraContenedor.aristas_internas.append(agujero[0])

            for arista in agujero:
                arista.cara = caraContenedor




