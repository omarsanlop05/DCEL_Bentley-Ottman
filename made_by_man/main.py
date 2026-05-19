"""
Algoritmo fusión de listas de aristas/fusión de caras.

Para hacer que interactuen las figuras, primero se debe correr el algoritmo de interseccion.
1. Buscar intersección
2. Se anotan las aristas que interactúan
3. Se remueven las aristas que interactúan de la lista de aristas
4. Por cada arista vamos a agregar dos aristas
	1. Los primos inician en el vertice original
	2. Y los primos primos, empiezan en X (El nuevo punto)
5. Hay que ordenarlo por ordenamiento circular
	1. Sacas los vectores con respecto a X, y los ordenas en círculo. Esto es un ordenamietno polar
	2. La pareja del primo es el primo primo del que era su pareja anterior, y la pareja del primo primo es el primo de su pareja anterior
	3. El siguiente de los primos primos no cambia, y el de los primos es el siguiente en el ordenamiento circular (el siguiente primo primo)
    4. El anterior de los primos sigue siendo el siguiente (en caso de que no exista, es el primo primo equivalente), y el de los primo primo es el anterior en la lista de ordenamiento circular (el anterior de los primos)
6. Rearme de caras

Una vez teniendo todas las intersecciones y reacomodando las aristas, el orden de los ciclos debería quedar automático.
Ahora se deben indentificar las caras, y para esto se utilizando los ciclos.

1. Extracción de ciclos (según el orden ya adquirido automáticamente.
2. Clasificar si son internos o externos.
    - Iniciando desde el vértice más a la izquierda dle ciclo:
        1. > 180° ley de la mano derecha (hacia la izquierda) --> Externo
        2. < 180° ley de la mano derecha (hacia la izquierda) --> Interno
3. Grafos de caras:
    -> ciclo = nodo
    -> se conectan los nodos si hay una arista a la izquierda del nodo más izquierdo
    -> sino, se conecta al infinito
4. Cada componente conectado se registra como cara
    - Se asigna como liga externa a una arista del ciclo interno
    - Se genera una lista de ligas internas con una arista de cada ciclo externo
5. Se actualiza cada arista con cara a la que pertenece

"""

"""
Algoritmo fusión de listas de aristas/fusión de caras.
"""

from Figuras import *
from ArbolBarrido import *
from Segmentacion import *
import heapq
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from PyGame import VisorDCEL

# ─────────────────────────────────────────────────────────────
# Estado global del barrido
# ─────────────────────────────────────────────────────────────
eventos = []
eventos_por_punto = {}
vistos = set()
intersecciones = {}
arbol = ArbolBarrido()


# ─────────────────────────────────────────────────────────────
# PASO 1 – Bentley-Ottmann
# ─────────────────────────────────────────────────────────────

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
        ev = Evento(inter, None, [sI, sD])
        heapq.heappush(eventos, ev)
        eventos_por_punto[inter] = ev
        intersecciones[clave] = {"punto": inter, "segmentos": [sI, sD]}
    else:
        if inter not in eventos_por_punto:
            return
        even = eventos_por_punto[inter]
        if sI not in even.contiene: even.agregarContiene(sI)
        if sD not in even.contiene: even.agregarContiene(sD)
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
        todos = list({id(s): s for s in L + U + C}.values())
        if clave not in intersecciones:
            intersecciones[clave] = {"punto": p, "segmentos": todos}
        else:
            for seg in todos:
                if seg not in intersecciones[clave]["segmentos"]:
                    intersecciones[clave]["segmentos"].append(seg)
        vistos.add(p)
    for seg in L + C:
        if arbol.buscar(seg): arbol.eliminar(seg)
    arbol.sweep_y = p.y - 1e-9
    horizontales = [s for s in U + C if abs(s.p1.y - s.p2.y) < 1e-9]
    no_horizontales = [s for s in U + C if abs(s.p1.y - s.p2.y) >= 1e-9]
    for seg in no_horizontales: arbol.insertar(seg)
    for seg in horizontales:
        _manejar_horizontal(seg)
        arbol.insertar(seg)
    if not U + C:
        segs = arbol.en_orden()
        indice, sI = 0, None
        while indice < len(segs) and arbol._x_en_sweep(segs[indice]) < p.x:
            sI = segs[indice];
            indice += 1
        sD = segs[indice] if indice < len(segs) else None
        encuentraEventos(sI, sD, p)
    else:
        UC = [s for s in arbol.en_orden() if s in U + C]
        if UC:
            for i in range(len(UC) - 1):
                encuentraEventos(UC[i], UC[i + 1], p)
            for seg in UC:
                izq, der = arbol.vecinos(seg)
                if izq and izq not in UC: encuentraEventos(izq, seg, p)
                if der and der not in UC: encuentraEventos(seg, der, p)


def encuentraIntersecciones(segmentos_lista):
    for segmento in segmentos_lista:
        if segmento.p1 not in vistos:
            vistos.add(segmento.p1)
            ev = Evento(segmento.p1, [segmento], None)
            heapq.heappush(eventos, ev)
            eventos_por_punto[segmento.p1] = ev
        else:
            eventos_por_punto[segmento.p1].agregarInicial(segmento)
        if segmento.p2 not in vistos:
            vistos.add(segmento.p2)
            ev = Evento(segmento.p2, None, None)
            heapq.heappush(eventos, ev)
            eventos_por_punto[segmento.p2] = ev
    while eventos:
        procesarEvento(heapq.heappop(eventos))
    return intersecciones


# ─────────────────────────────────────────────────────────────
# NUEVO PASO 2, 3 y 4: RECONSTRUCCIÓN LIMPIA Y ALGORITMO DE CARAS
# ─────────────────────────────────────────────────────────────

def area_con_signo(ciclo):
    area = 0
    for a in ciclo:
        p1 = a.origen.pt
        p2 = a.antiarista.origen.pt
        area += (p1.x * p2.y) - (p2.x * p1.y)
    return area / 2.0


def punto_en_poligono(pt, ciclo):
    cruces = 0
    for a in ciclo:
        p1 = a.origen.pt
        p2 = a.antiarista.origen.pt
        if min(p1.y, p2.y) <= pt.y < max(p1.y, p2.y):
            ix = p1.x + (pt.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y)
            if ix > pt.x:
                cruces += 1
    return cruces % 2 == 1


def _punto_en_segmento(pt, seg, tol=1e-2):
    import math
    dist_num = abs(
        (seg.p2.y - seg.p1.y) * pt.x - (seg.p2.x - seg.p1.x) * pt.y + seg.p2.x * seg.p1.y - seg.p2.y * seg.p1.x)
    longitud = math.hypot(seg.p2.x - seg.p1.x, seg.p2.y - seg.p1.y)
    if longitud == 0: return False

    if (dist_num / longitud) < tol:
        min_x, max_x = min(seg.p1.x, seg.p2.x), max(seg.p1.x, seg.p2.x)
        min_y, max_y = min(seg.p1.y, seg.p2.y), max(seg.p1.y, seg.p2.y)
        if (min_x - tol <= pt.x <= max_x + tol) and (min_y - tol <= pt.y <= max_y + tol):
            return True
    return False


# MODIFICACIÓN CLAVE: Red de seguridad universal
def preprocesar_geometria(segmentos_globales, tolerancia_decimales=2):
    """
    Actúa como red de seguridad obligando a revisar cruces ignorados,
    garantizando que capas como el ARBUSTO resuelvan sus propias auto-intersecciones.
    """
    import itertools
    from Figuras import Punto

    # 1. Snap to grid (Alineación de geometría)
    for seg in segmentos_globales:
        seg.p1.x = round(seg.p1.x, tolerancia_decimales)
        seg.p1.y = round(seg.p1.y, tolerancia_decimales)
        seg.p2.x = round(seg.p2.x, tolerancia_decimales)
        seg.p2.y = round(seg.p2.y, tolerancia_decimales)

        for inter in seg.intersecciones:
            inter.x = round(inter.x, tolerancia_decimales)
            inter.y = round(inter.y, tolerancia_decimales)

    # 2. Búsqueda Exhaustiva Total (Fuerza bruta para no omitir nada)
    for sA, sB in itertools.combinations(segmentos_globales, 2):

        # A. Atrapa cruces "X" que Bentley-Ottmann haya perdido
        pts = sA.interseccion(sB)
        if pts:
            for p_inter in pts:
                p_red = Punto(round(p_inter.x, tolerancia_decimales), round(p_inter.y, tolerancia_decimales))
                if not any(p.comparar(p_red) for p in sA.intersecciones):
                    sA.intersecciones.append(p_red)
                if not any(p.comparar(p_red) for p in sB.intersecciones):
                    sB.intersecciones.append(p_red)

        # B. Atrapa "T-Junctions" y colineales (Punta contra Cuerpo)
        for pt in [sB.p1, sB.p2]:
            if not pt.comparar(sA.p1) and not pt.comparar(sA.p2):
                if _punto_en_segmento(pt, sA):
                    nuevo_pt = Punto(pt.x, pt.y)
                    if not any(p.comparar(nuevo_pt) for p in sA.intersecciones):
                        sA.intersecciones.append(nuevo_pt)

        for pt in [sA.p1, sA.p2]:
            if not pt.comparar(sB.p1) and not pt.comparar(sB.p2):
                if _punto_en_segmento(pt, sB):
                    nuevo_pt = Punto(pt.x, pt.y)
                    if not any(p.comparar(nuevo_pt) for p in sB.intersecciones):
                        sB.intersecciones.append(nuevo_pt)


def reconstruir_overlay_y_caras(segmentos_globales):
    print("\n── Paso 2 y 3: Reconstrucción Limpia de la Topología DCEL ──")

    # 1. DE-DUPLICACIÓN MATEMÁTICA
    sub_segmentos_unicos = set()
    for seg in segmentos_globales:
        pts = [seg.p1, seg.p2] + seg.intersecciones
        pts.sort(key=lambda p: (round(p.x, 6), round(p.y, 6)))

        for i in range(len(pts) - 1):
            pA, pB = pts[i], pts[i + 1]
            tA = (round(pA.x, 6), round(pA.y, 6))
            tB = (round(pB.x, 6), round(pB.y, 6))
            if tA != tB:
                sub_segmentos_unicos.add((tA, tB) if tA < tB else (tB, tA))

    print(f"  > Sub-segmentos únicos filtrados: {len(sub_segmentos_unicos)}")

    # 2. CREACIÓN DCEL BÁSICA
    vertices_g = {}
    aristas_g = {}

    for tA, tB in sub_segmentos_unicos:
        if tA not in vertices_g: vertices_g[tA] = Vertice(f"V_{len(vertices_g)}", Punto(tA[0], tA[1]))
        if tB not in vertices_g: vertices_g[tB] = Vertice(f"V_{len(vertices_g)}", Punto(tB[0], tB[1]))

        vA, vB = vertices_g[tA], vertices_g[tB]

        a1 = Arista(f"E_{len(aristas_g)}")
        a2 = Arista(f"E_{len(aristas_g) + 1}")

        a1.origen, a2.origen = vA, vB
        a1.antiarista, a2.antiarista = a2, a1

        aristas_g[a1.nombre] = a1
        aristas_g[a2.nombre] = a2

        if not hasattr(vA, 'outgoing'): vA.outgoing = []
        if not hasattr(vB, 'outgoing'): vB.outgoing = []

        vA.outgoing.append(a1)
        vB.outgoing.append(a2)
        vA.arista_adyacente = a1
        vB.arista_adyacente = a2

    # 3. GRAFO RADIAL: Enlazar Next y Prev
    for v in vertices_g.values():
        v.outgoing.sort(key=lambda a: math.atan2(a.antiarista.origen.pt.y - v.pt.y,
                                                 a.antiarista.origen.pt.x - v.pt.x))
        k = len(v.outgoing)
        for i in range(k):
            e_out = v.outgoing[i]
            e_in = e_out.antiarista
            e_next = v.outgoing[(i - 1) % k]

            e_in.siguiente = e_next
            e_next.anterior = e_in

    print(f"  > DCEL topológicamente perfecta: {len(vertices_g)} vértices, {len(aristas_g)} medias-aristas")

    print("\n── Paso 4: Rearme de Caras y Jerarquías ──")

    # 4. EXTRACCIÓN DE CICLOS
    ciclos = []
    visitadas = set()
    for arista in aristas_g.values():
        if arista.nombre not in visitadas:
            ciclo = []
            actual = arista
            while actual.nombre not in visitadas:
                visitadas.add(actual.nombre)
                ciclo.append(actual)
                actual = actual.siguiente
            ciclos.append(ciclo)

    # 5. CLASIFICACIÓN
    exteriores = []
    agujeros = []
    for ciclo in ciclos:
        area = area_con_signo(ciclo)
        if area > 0:
            exteriores.append((ciclo, area))
        else:
            agujeros.append((ciclo, area))

    print(f"  > Ciclos hallados: {len(exteriores)} Fronteras Externas, {len(agujeros)} Huecos")

    # 6. INSTANCIACIÓN DE CARAS Y ASIGNACIÓN DE HUECOS
    caras_g = {}
    cara_infinita = Cara("f_infinita")
    caras_g["f_infinita"] = cara_infinita

    for i, (ext, _) in enumerate(exteriores):
        f = Cara(f"f{i}")
        f.aristas_externas = ext[0]
        f.aristas_internas = []
        for a in ext:
            a.cara = f
        caras_g[f.nombre] = f
        ext[0]._cara_obj = f

    for agujero, _ in agujeros:
        v_izq = min(agujero, key=lambda a: a.origen.pt.x).origen.pt
        padres_candidatos = []

        for ext, area in exteriores:
            if punto_en_poligono(v_izq, ext):
                padres_candidatos.append((ext, area))

        if padres_candidatos:
            padre_directo = min(padres_candidatos, key=lambda x: x[1])[0]
            f = padre_directo[0]._cara_obj
            f.aristas_internas.append(agujero[0])
            for a in agujero: a.cara = f
        else:
            cara_infinita.aristas_internas.append(agujero[0])
            for a in agujero: a.cara = cara_infinita

    print(f"  > Total final: {len(caras_g)} caras válidas ensambladas (incluyendo f_infinita)")
    return vertices_g, aristas_g, caras_g


# ─────────────────────────────────────────────────────────────
# EJECUCIÓN PRINCIPAL
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Asegúrate de colocar las rutas correctas a los archivos
    archivos = ["./Proyecto_SubdivisionesInteractivas/Arbusto/layerARBUSTO" + str(i) for i in range(1, 3)]

    print("── Paso 1: Intersecciones (Bentley-Ottmann) ──")
    segmentos_globales = extraer_todos_los_segmentos(archivos)
    diccionario_intersecciones = encuentraIntersecciones(segmentos_globales)
    print(
        f"  > {len(segmentos_globales)} segmentos base, {len(diccionario_intersecciones)} ptos de cruce hallados por Sweep-Line")

    for clave, info in diccionario_intersecciones.items():
        punto = info["punto"]
        for seg in info["segmentos"]:
            if not punto.comparar(seg.p1) and not punto.comparar(seg.p2):
                if punto not in seg.intersecciones:
                    seg.intersecciones.append(punto)

    # AQUÍ ESTÁ LA MAGIA: Obligamos al código a repasar TODO por si acaso
    print("── Pre-procesado: Red de seguridad para cruces internos y T-Junctions ──")
    preprocesar_geometria(segmentos_globales, tolerancia_decimales=2)

    # Correr nuestro motor de reconstrucción consolidado
    vertices_g, aristas_g, caras_g = reconstruir_overlay_y_caras(segmentos_globales)

    ##Archivos de vertices, aristas y caras
    #######################################
    with open("vertices.txt", "w", encoding="utf-8") as f:
        f.write("Archivo de vértices\n")
        f.write("#" * 30 + "\n")
        f.write(f"{'Nombre':<10} {'x':<10} {'y':<10} {'Incidente':<10}\n")
        f.write("#" * 30 + "\n")

        for _, v in vertices_g.items():
            nombre = v.nombre
            x = v.pt.x
            y = v.pt.y
            incidente = v.arista_adyacente.nombre if v.arista_adyacente else "None"

            f.write(f"{nombre:<10} {x:<10} {y:<10} {incidente:<10}\n")

    with open("aristas.txt", "w", encoding="utf-8") as f:
        # Encabezado
        f.write("Archivo de aristas\n")
        f.write("#" * 40 + "\n")
        f.write(f"{'Nombre':<10} {'Origen':<10} {'Pareja':<10} {'Cara':<10} {'Sigue':<10} {'Antes':<10}\n")
        f.write("#" * 40 + "\n")

        for nombre, a in aristas_g.items():
            origen = a.origen.nombre if a.origen else "None"
            pareja = a.antiarista.nombre if a.antiarista else "None"
            cara = a.cara.nombre if a.cara else "None"
            siguiente = a.siguiente.nombre if a.siguiente else "None"
            anterior = a.anterior.nombre if a.anterior else "None"

            f.write(f"{nombre:<10} {origen:<10} {pareja:<10} {cara:<10} {siguiente:<10} {anterior:<10}\n")

    with open("caras.txt", "w", encoding="utf-8") as f:
        f.write("Archivo de caras\n")
        f.write("#" * 22 + "\n")
        f.write(f"{'Nombre':<10} {'Interno':<12} {'Externo':<10}\n")
        f.write("#" * 22 + "\n")

        for nombre, c in caras_g.items():
            externo = c.aristas_externas.nombre if c.aristas_externas else "None"

            # Las internas pueden ser varias, las ponemos entre corchetes
            if c.aristas_internas:
                interno = "[" + ", ".join(a.nombre for a in c.aristas_internas) + "]"
            else:
                interno = "None"

            f.write(f"{nombre:<10} {interno:<12} {externo:<10}\n")


    # ── Visualización Mejorada ──
    fig, ax = plt.subplots(figsize=(10, 10))
    colores_caras = cm.Set3.colors

    for i, (nombre, cara) in enumerate(caras_g.items()):
        if nombre == "f_infinita" or cara.aristas_externas is None:
            continue

        puntos = []
        actual = cara.aristas_externas
        visitadas = set()
        while actual and actual.nombre not in visitadas:
            visitadas.add(actual.nombre)
            puntos.append((actual.origen.pt.x, actual.origen.pt.y))
            actual = actual.siguiente

        if len(puntos) >= 3:
            xs = [p[0] for p in puntos]
            ys = [p[1] for p in puntos]
            color = colores_caras[i % len(colores_caras)] + (0.6,)
            ax.fill(xs, ys, color=color, zorder=1)
            cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
            ax.text(cx, cy, nombre, fontsize=9, ha='center', va='center',
                    fontweight='bold', zorder=4, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    dibujadas = set()
    for nombre, a in aristas_g.items():
        if nombre in dibujadas: continue
        p1 = a.origen.pt
        p2 = a.antiarista.origen.pt
        ax.plot([p1.x, p2.x], [p1.y, p2.y], 'k-', linewidth=1.5, alpha=0.9, zorder=3)
        dibujadas.add(nombre)
        dibujadas.add(a.antiarista.nombre)

    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.4)
    plt.title(f"Resultado Final DCEL: {len(caras_g)} Caras Identificadas")
    plt.tight_layout()
    plt.show()


VisorDCEL(vertices_g, aristas_g, caras_g).ejecutar()

