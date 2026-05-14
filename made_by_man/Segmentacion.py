import Figuras

# Heredamos de Figuras.Segmento para que el algoritmo de barrido
# (Bentley-Ottmann) lo procese sin problemas y pueda usar sus métodos matemáticos.
class SegmentoGeometrico(Figuras.Segmento):
    def __init__(self, p1, p2, layer_id):
        super().__init__(p1, p2) # Esto acomoda p1 y p2 con la lógica de Figuras.py
        self.layer_id = layer_id
        # Aquí guardaremos los puntos donde este segmento choca con otros
        self.intersecciones = []

    def __str__(self):
        return f"Layer {self.layer_id} | ({self.p1.x}, {self.p1.y}) -> ({self.p2.x}, {self.p2.y})"


def extraer_todos_los_segmentos(rutas_archivos):
    segmentos_totales = []

    for layer_id, ruta in enumerate(rutas_archivos):
        # Usamos tu función definirObjetos
        vertices, aristas, caras = Figuras.definirObjetos(ruta)

        visitadas = set()

        for a in aristas.values():
            if a.nombre in visitadas:
                continue

            visitadas.add(a.nombre)
            if a.antiarista:
                visitadas.add(a.antiarista.nombre)

            # En tu Vertice, el punto se llama 'pt', no 'punto'
            p1 = a.origen.pt
            p2 = a.antiarista.origen.pt

            # Guardamos el segmento independiente de la DCEL
            nuevo_segmento = SegmentoGeometrico(p1, p2, layer_id)
            segmentos_totales.append(nuevo_segmento)

    return segmentos_totales