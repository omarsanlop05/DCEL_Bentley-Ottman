

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

