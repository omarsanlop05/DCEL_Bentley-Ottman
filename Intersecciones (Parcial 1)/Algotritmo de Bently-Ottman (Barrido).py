from Punto import Punto
from Segmento import Segmento
from Evento import Evento
from ArbolBarrido import ArbolBarrido
import matplotlib.pyplot as plt
from matplotlib import animation
import heapq
import time

inicio = time.perf_counter()

segmentos = []
nombres_segmentos = {}

eventos = []
eventos_por_punto = {}

vistos = set()

intersecciones = {}

arbol = ArbolBarrido()

X1 = []
Y1 = []
X2 = []
Y2 = []

with open("0.txt") as f:
    for linea in f:
        x1, y1, x2, y2, nombre = linea.strip().split()

        X1.append(int(x1))
        X2.append(int(x2))
        Y1.append(int(y1))
        Y2.append(int(y2))

        p1 = Punto(int(x1), int(y1))
        p2 = Punto(int(x2), int(y2))

        if p1.y < p2.y or (p1.y == p2.y and p1.x > p2.x):
            p1, p2 = p2, p1

        s = Segmento(p1, p2)
        num_segmentos = len(segmentos)

        segmentos.append(s)
        nombres_segmentos[id(s)] = f"S{num_segmentos}"

def encuentraEventos(sI, sD, p):
    if sI is None or sD is None:
        return

    pts = sI.interseccion(sD)
    if not pts:
        return

    inter = pts[0]
    clave = (round(inter.x, 9), round(inter.y, 9))

    debajo = inter.y < p.y - 1e-9
    mismo_nivel_der = (abs(inter.y - p.y) < 1e-9 and inter.x >= p.x - 1e-9)

    if not (debajo or mismo_nivel_der):
        return

    if inter not in vistos:
        vistos.add(inter)
        eventoInter = Evento(inter, None, [sI, sD])
        heapq.heappush(eventos, eventoInter)
        eventos_por_punto[inter] = eventoInter
        intersecciones[clave] = {"punto": inter, "segmentos": [sI, sD]}
    else:
        if inter not in eventos_por_punto:
            return
        even = eventos_por_punto[inter]
        if sI not in even.contiene:
            even.agregarContiene(sI)
        if sD not in even.contiene:
            even.agregarContiene(sD)

        if clave in intersecciones:
            for seg in [sI, sD]:
                if seg not in intersecciones[clave]["segmentos"]:
                    intersecciones[clave]["segmentos"].append(seg)

def _manejar_horizontal(seg_h):
    x_min = min(seg_h.p1.x, seg_h.p2.x)
    x_max = max(seg_h.p1.x, seg_h.p2.x)

    for seg in arbol.en_orden():
        x_seg = arbol._x_en_sweep(seg)
        if x_min - 1e-9 <= x_seg <= x_max + 1e-9:
            pts = seg_h.interseccion(seg)
            for inter in pts:
                clave = (round(inter.x, 9), round(inter.y, 9))
                if clave not in intersecciones:
                    intersecciones[clave] = {"punto": inter, "segmentos": [seg_h, seg]}
                else:
                    for s in [seg_h, seg]:
                        if s not in intersecciones[clave]["segmentos"]:
                            intersecciones[clave]["segmentos"].append(s)

def procesarEvento(evento):
    p = evento.punto
    U = evento.inicio
    C = evento.contiene
    L = [seg for seg in arbol.en_orden() if seg.p2 == p]

    C = list(set(C) - set(U) - set(L))

    if len(U + C + L) > 1:
        clave = (round(p.x, 9), round(p.y, 9))
        todos = list({id(s): s for s in L + U + C}.values())  # deduplicar
        if clave not in intersecciones:
            intersecciones[clave] = {"punto": p, "segmentos": todos}
        else:
            for seg in todos:
                if seg not in intersecciones[clave]["segmentos"]:
                    intersecciones[clave]["segmentos"].append(seg)
        vistos.add(p)

    for seg in L + C:
        if arbol.buscar(seg):
            arbol.eliminar(seg)

    arbol.sweep_y = p.y - 1e-9

    horizontales = [seg for seg in U + C if abs(seg.p1.y - seg.p2.y) < 1e-9]
    no_horizontales = [seg for seg in U + C if abs(seg.p1.y - seg.p2.y) >= 1e-9]

    for seg in no_horizontales:
        arbol.insertar(seg)

    for seg in horizontales:
        _manejar_horizontal(seg)
        arbol.insertar(seg)

    if len(U+C) == 0:
        indice = 0
        segs = arbol.en_orden()
        sI = None

        while indice < len(segs) and arbol._x_en_sweep(segs[indice]) < p.x:
            sI = segs[indice]
            indice+=1

        sD = None
        if indice < len(segs):
            sD = segs[indice]

        encuentraEventos(sI, sD, p)

    else:
        UC = [seg for seg in arbol.en_orden() if seg in U + C]
        if UC:
            for i in range(len(UC) - 1):
                encuentraEventos(UC[i], UC[i + 1], p)

            for seg in UC:
                izq, der = arbol.vecinos(seg)
                if izq and izq not in UC:
                    encuentraEventos(izq, seg, p)
                if der and der not in UC:
                    encuentraEventos(seg, der, p)

def encuentraIntersecciones(segmentos):
    for segmento in segmentos:

        if segmento.p1 not in vistos:
            vistos.add(segmento.p1)
            eventoInicial = Evento(segmento.p1, [segmento], None)
            heapq.heappush(eventos, eventoInicial)
            eventos_por_punto[segmento.p1] = eventoInicial
        else:
            even = eventos_por_punto[segmento.p1]
            even.agregarInicial(segmento)

        if segmento.p2 not in vistos:
            vistos.add(segmento.p2)
            eventoFinal = Evento(segmento.p2, None, None)
            heapq.heappush(eventos, eventoFinal)
            eventos_por_punto[segmento.p2] = eventoFinal

    while eventos:
        evento = heapq.heappop(eventos)
        procesarEvento(evento)

    return intersecciones

ints = encuentraIntersecciones(segmentos)
print("Intersecciones encontradas: ", len(ints))

fin = time.perf_counter()
print(f"Tiempo de ejecución: {fin - inicio:.6f} segundos")

xI = []
yI = []

for inter in ints.values():
    p = inter["punto"]
    xI.append(p.x)
    yI.append(p.y)

fig = plt.figure()
ax = fig.add_subplot()

xmax = max(X1 + X2)
xmin = min(X1 + X2)
ymax = max(Y1 + Y2)
ymin = min(Y1 + Y2)

fig, ax = plt.subplots()
ax.set_xlim(xmin - 2, xmax + 2)
ax.set_ylim(ymin - 2, ymax + 2)

for spine in ax.spines:
    ax.spines[spine].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

lineas = []
for s in segmentos:
    linea_seg, = ax.plot((s.p1.x, s.p2.x), (s.p1.y, s.p2.y), color="black", zorder=-1)
    lineas.append((s, linea_seg))

ax.scatter([s.p1.x for s in segmentos] + [s.p2.x for s in segmentos],
           [s.p1.y for s in segmentos] + [s.p2.y for s in segmentos],
           color="red", zorder=2)

if xI:
    ax.scatter(xI, yI, color="yellow", zorder=3)

sweep_line = ax.axhline(y=ymax, color="green", linewidth=2)


def cruza(segmento, y0):
    y1, y2 = segmento.p1.y, segmento.p2.y
    return (y1 - y0) * (y2 - y0) <= 0


def anim(frame):
    linea_y = ymax - frame * (ymax - ymin) / n

    for s, linea_seg in lineas:
        if cruza(s, linea_y):
            linea_seg.set_color("red")
        else:
            linea_seg.set_color("black")

    sweep_line.set_ydata([linea_y, linea_y])
    return [sweep_line] + [l for _, l in lineas]


print("Creando animación")
n = 120
animacion = animation.FuncAnimation(fig, anim, frames=n, blit=True)
html_content = animacion.to_jshtml()

with open("animBarrido.html", "w", encoding="utf-8") as file:
    file.write(html_content)

print("Animacion creada correctamente")

fin = time.perf_counter()
print(f"Tiempo de ejecución: {fin - inicio:.6f} segundos")