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
Algoritmo fusión de listas de aristas / fusión de caras.

ESTADO ACTUAL → Paso 2:
  ✅ Paso 1 – Buscar intersecciones (Bentley-Ottmann ya funcionando)
  ✅ Paso 2 – Cargar DCEL unificada (con prefijos por capa para evitar
              colisiones de nombres) y mapear cada SegmentoGeométrico
              a su(s) semiarista(s) DCEL correspondiente(s)
  🔜 Paso 3 – Subdivisión: primas y primas-primas por cada intersección
  🔜 Paso 4 – Ordenamiento polar + reasignación de siguiente/anterior/pareja
  🔜 Paso 5 – Rearme de caras
"""

from Figuras import *
from Graficacion import *
from ArbolBarrido import *
from Segmentacion import *
import heapq
import math
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────
# Estado global del barrido
# ─────────────────────────────────────────────────────────────
eventos            = []
eventos_por_punto  = {}
vistos             = set()
intersecciones     = {}
arbol              = ArbolBarrido()

# ─────────────────────────────────────────────────────────────
# Bentley-Ottmann (sin cambios)
# ─────────────────────────────────────────────────────────────

def encuentraEventos(sI, sD, p):
    if sI is None or sD is None:
        return
    pts = sI.interseccion(sD)
    if not pts:
        return
    inter = pts[0]
    clave = (round(inter.x, 9), round(inter.y, 9))
    debajo          = inter.y < p.y - 1e-9
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
    horizontales    = [s for s in U + C if abs(s.p1.y - s.p2.y) < 1e-9]
    no_horizontales = [s for s in U + C if abs(s.p1.y - s.p2.y) >= 1e-9]
    for seg in no_horizontales: arbol.insertar(seg)
    for seg in horizontales:
        _manejar_horizontal(seg)
        arbol.insertar(seg)
    if not U + C:
        segs = arbol.en_orden()
        indice, sI = 0, None
        while indice < len(segs) and arbol._x_en_sweep(segs[indice]) < p.x:
            sI = segs[indice]; indice += 1
        sD = segs[indice] if indice < len(segs) else None
        encuentraEventos(sI, sD, p)
    else:
        UC = [s for s in arbol.en_orden() if s in U + C]
        if UC:
            for i in range(len(UC) - 1):
                encuentraEventos(UC[i], UC[i+1], p)
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
# PASO 2 – Cargar DCEL unificada con prefijos por capa
# ─────────────────────────────────────────────────────────────

def cargar_dcel_unificada(archivos):
    """
    Carga la DCEL de cada capa y las combina en diccionarios globales.

    PROBLEMA RESUELTO: capas distintas pueden tener vértices/aristas con el
    mismo nombre (p.ej. 'p1' en layer01 y 'p1' en layer03 son puntos distintos).
    Solución: prefijar cada nombre con 'L{i}_' según el índice de la capa.

    El renombrado se propaga a todas las referencias internas de la DCEL
    (origen, antiarista, cara, siguiente, anterior, arista_adyacente).
    """
    vertices_g = {}
    aristas_g  = {}
    caras_g    = {}

    for layer_id, ruta in enumerate(archivos):
        v, a, c = definirObjetos(ruta)
        pre = f"L{layer_id}_"

        # 1. Renombrar nombres propios
        for obj in list(v.values()) + list(a.values()) + list(c.values()):
            obj.nombre = pre + obj.nombre

        # 2. Insertar en dicts globales con nueva clave
        for nombre, vertice in v.items():
            vertices_g[pre + nombre] = vertice
        for nombre, arista in a.items():
            aristas_g[pre + nombre] = arista
        for nombre, cara in c.items():
            caras_g[pre + nombre] = cara

        # 3. Actualizar referencias internas (los objetos ya tienen nombre prefijado,
        #    así que buscamos directamente en los dicts globales)
        for arista in a.values():
            if arista.origen:
                arista.origen     = vertices_g.get(arista.origen.nombre)
            if arista.antiarista:
                arista.antiarista = aristas_g.get(arista.antiarista.nombre)
            if arista.cara:
                arista.cara       = caras_g.get(arista.cara.nombre)
            if arista.siguiente:
                arista.siguiente  = aristas_g.get(arista.siguiente.nombre)
            if arista.anterior:
                arista.anterior   = aristas_g.get(arista.anterior.nombre)

        for vertice in v.values():
            if vertice.arista_adyacente:
                vertice.arista_adyacente = aristas_g.get(vertice.arista_adyacente.nombre)

        for cara in c.values():
            cara.aristas_internas = [
                aristas_g.get(ar.nombre, ar) for ar in cara.aristas_internas
            ]
            if cara.aristas_externas:
                cara.aristas_externas = aristas_g.get(cara.aristas_externas.nombre)

    print(f"  DCEL unificada: {len(vertices_g)} vértices, "
          f"{len(aristas_g)} aristas, {len(caras_g)} caras")
    return vertices_g, aristas_g, caras_g


def mapear_segmentos_a_aristas(segmentos_geo, aristas_g):
    """
    Para cada SegmentoGeométrico encuentra las semiaristas DCEL que
    representan geométricamente ese mismo segmento (en cualquier dirección).

    Retorna:  id(segmento) → [Arista, ...]
    """
    def pts_iguales(p1, p2):
        return math.isclose(p1.x, p2.x, abs_tol=1e-7) and \
               math.isclose(p1.y, p2.y, abs_tol=1e-7)

    mapa = {}
    for seg in segmentos_geo:
        encontradas = []
        for arista in aristas_g.values():
            orig = arista.origen.pt
            dest = arista.antiarista.origen.pt
            if ((pts_iguales(orig, seg.p1) and pts_iguales(dest, seg.p2)) or
                (pts_iguales(orig, seg.p2) and pts_iguales(dest, seg.p1))):
                encontradas.append(arista)
        mapa[id(seg)] = encontradas

    sin_match = [s for s in segmentos_geo if not mapa[id(s)]]
    print(f"  Segmentos mapeados a DCEL: "
          f"{len(segmentos_geo) - len(sin_match)}/{len(segmentos_geo)}")
    if sin_match:
        print(f"  ⚠️  Sin match ({len(sin_match)}):")
        for s in sin_match:
            print(f"      Layer {s.layer_id}: {s.p1} → {s.p2}")
    return mapa


def filtrar_intersecciones_entre_capas(diccionario_intersecciones):
    """
    Conserva solo intersecciones entre segmentos de capas DISTINTAS.
    Los vértices compartidos dentro de una misma capa ya están en la DCEL.
    """
    reales = {
        clave: info
        for clave, info in diccionario_intersecciones.items()
        if len(set(s.layer_id for s in info["segmentos"])) > 1
    }
    print(f"  Intersecciones entre capas: {len(reales)} "
          f"(de {len(diccionario_intersecciones)} totales)")
    return reales


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":

    archivos = ["./layers/layer0" + str(i) for i in range(1, 6)]

    # ── Paso 1: Barrido geométrico ───────────────────────────
    print("── Paso 1: Extracción de segmentos e intersecciones ──")
    segmentos_globales = extraer_todos_los_segmentos(archivos)
    print(f"  Segmentos extraídos: {len(segmentos_globales)}")

    diccionario_intersecciones = encuentraIntersecciones(segmentos_globales)
    print(f"  Puntos de intersección encontrados: {len(diccionario_intersecciones)}")

    for clave, info in diccionario_intersecciones.items():
        punto = info["punto"]
        for seg in info["segmentos"]:
            if not punto.comparar(seg.p1) and not punto.comparar(seg.p2):
                if punto not in seg.intersecciones:
                    seg.intersecciones.append(punto)

    # ── Paso 2: DCEL unificada ───────────────────────────────
    print("\n── Paso 2: Carga y unificación de la DCEL ──")
    vertices_g, aristas_g, caras_g = cargar_dcel_unificada(archivos)
    mapa_seg_arista                 = mapear_segmentos_a_aristas(segmentos_globales, aristas_g)
    intersecciones_reales           = filtrar_intersecciones_entre_capas(diccionario_intersecciones)

    # Verificación: mostrar aristas DCEL involucradas
    aristas_a_subdividir = set()
    for info in intersecciones_reales.values():
        for seg in info["segmentos"]:
            for arista in mapa_seg_arista.get(id(seg), []):
                aristas_a_subdividir.add(arista.nombre)

    print(f"\n  Semiaristas DCEL a subdividir: {len(aristas_a_subdividir)}")

    # Verificación de integridad de la DCEL unificada
    print("\n  Verificación de integridad:")
    errores = 0
    for nombre, a in aristas_g.items():
        if a.origen is None:       print(f"    ⚠️  {nombre}: origen None");       errores += 1
        if a.antiarista is None:   print(f"    ⚠️  {nombre}: antiarista None");   errores += 1
        if a.siguiente is None:    print(f"    ⚠️  {nombre}: siguiente None");    errores += 1
        if a.anterior is None:     print(f"    ⚠️  {nombre}: anterior None");     errores += 1
    if errores == 0:
        print("    ✅ Todas las referencias de la DCEL son válidas")
    else:
        print(f"    ❌ {errores} referencias inválidas")

    # ── Visualización ────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 10))
    colores = ['blue', 'green', 'orange', 'purple', 'brown']

    for seg in segmentos_globales:
        color = colores[seg.layer_id % len(colores)]
        ax.plot([seg.p1.x, seg.p2.x], [seg.p1.y, seg.p2.y],
                color=color, linewidth=1.5, alpha=0.7,
                label=f"Layer {seg.layer_id + 1}")
        for pt in seg.intersecciones:
            ax.plot(pt.x, pt.y, 'ro', markersize=5, zorder=5)

    for v in vertices_g.values():
        ax.plot(v.pt.x, v.pt.y, 'ks', markersize=5, zorder=6)
        ax.text(v.pt.x + 0.1, v.pt.y + 0.1, v.nombre, fontsize=6, color='#333')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')

    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.4)
    plt.title("Paso 2 – DCEL unificada + intersecciones entre capas (puntos rojos)")
    plt.tight_layout()
    plt.savefig("paso2_dcel_unificada.png", dpi=150)
    plt.show()
