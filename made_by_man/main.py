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
from made_by_man.Caras import ConstructorCaras

"""
Algoritmo fusión de listas de aristas / fusión de caras.

ESTADO ACTUAL → Paso 3 (corregido):
  ✅ Paso 1 – Buscar intersecciones (Bentley-Ottmann)
  ✅ Paso 2 – DCEL unificada con prefijos por capa + mapeo segmento→arista
  ✅ Paso 3 – Subdivisión con DOS PASADAS:
              Pasada 1: crear X, primas y pp, registrar mapa global de
                        reemplazos (arista_eliminada → su_primo).
              Pasada 2: asignar siguiente/anterior usando el mapa global,
                        para que referencias a aristas ya eliminadas en
                        otra intersección se resuelvan correctamente.
  🔜 Paso 4 – Rearme de caras
"""

from Figuras import *
from Graficacion import *
from ArbolBarrido import *
from Segmentacion import *
import heapq
import math
import matplotlib.pyplot as plt

import matplotlib.cm as cm

# ─────────────────────────────────────────────────────────────
# Estado global del barrido
# ─────────────────────────────────────────────────────────────
eventos           = []
eventos_por_punto = {}
vistos            = set()
intersecciones    = {}
arbol             = ArbolBarrido()

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
# PASO 2 – DCEL unificada con prefijos por capa
# ─────────────────────────────────────────────────────────────

def cargar_dcel_unificada(archivos):
    """
    Carga la DCEL de cada capa y las combina en diccionarios globales.
    Prefixa cada nombre con 'L{i}_' para evitar colisiones entre capas.
    """
    vertices_g, aristas_g, caras_g = {}, {}, {}

    for layer_id, ruta in enumerate(archivos):
        v, a, c = definirObjetos(ruta)
        pre = f"L{layer_id}_"

        for obj in list(v.values()) + list(a.values()) + list(c.values()):
            obj.nombre = pre + obj.nombre

        for nombre, vertice in v.items(): vertices_g[pre + nombre] = vertice
        for nombre, arista  in a.items(): aristas_g [pre + nombre] = arista
        for nombre, cara    in c.items(): caras_g   [pre + nombre] = cara

        for arista in a.values():
            if arista.origen:     arista.origen     = vertices_g.get(arista.origen.nombre)
            if arista.antiarista: arista.antiarista = aristas_g .get(arista.antiarista.nombre)
            if arista.cara:       arista.cara       = caras_g   .get(arista.cara.nombre)
            if arista.siguiente:  arista.siguiente  = aristas_g .get(arista.siguiente.nombre)
            if arista.anterior:   arista.anterior   = aristas_g .get(arista.anterior.nombre)

        for vertice in v.values():
            if vertice.arista_adyacente:
                vertice.arista_adyacente = aristas_g.get(vertice.arista_adyacente.nombre)

        for cara in c.values():
            cara.aristas_internas = [aristas_g.get(ar.nombre, ar) for ar in cara.aristas_internas]
            if cara.aristas_externas:
                cara.aristas_externas = aristas_g.get(cara.aristas_externas.nombre)

    print(f"  {len(vertices_g)} vértices, {len(aristas_g)} aristas, {len(caras_g)} caras")
    return vertices_g, aristas_g, caras_g


def mapear_segmentos_a_aristas(segmentos_geo, aristas_g):
    """Conecta cada SegmentoGeométrico con sus semiaristas DCEL."""
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
    print(f"  Segmentos mapeados: {len(segmentos_geo)-len(sin_match)}/{len(segmentos_geo)}")
    if sin_match:
        for s in sin_match:
            print(f"    ⚠️  Layer {s.layer_id}: {s.p1}→{s.p2} sin match")
    return mapa


def filtrar_intersecciones_entre_capas(diccionario_intersecciones):
    """Conserva solo los puntos de intersección entre capas distintas."""
    reales = {
        k: v for k, v in diccionario_intersecciones.items()
        if len(set(s.layer_id for s in v["segmentos"])) > 1
    }
    print(f"  Intersecciones entre capas: {len(reales)} "
          f"(de {len(diccionario_intersecciones)} totales)")
    return reales

# ─────────────────────────────────────────────────────────────
# PASO 3 – Subdivisión en dos pasadas
# ─────────────────────────────────────────────────────────────

_contador_x = [0]
_contador_a = [0]


def _resolver(arista, mapa_eliminadas_a_primo, mapa_eliminadas_a_pp):
    """
    Dado un objeto Arista (posiblemente eliminado), devuelve la arista
    activa que lo reemplaza, o la misma si no fue eliminada.
    Útil para resolver referencias que quedaron apuntando a eliminadas.
    """
    nombre = arista.nombre if arista else None
    if nombre in mapa_eliminadas_a_primo:
        return mapa_eliminadas_a_primo[nombre]   # reemplazado por su primo
    return arista


def subdividir_todas(intersecciones_reales, mapa_seg_arista, vertices_g, aristas_g):
    """
    Subdivide todas las intersecciones en DOS PASADAS para resolver
    correctamente las referencias cruzadas entre intersecciones distintas.

    PASADA 1: por cada intersección, crear X + primas + pp.
              Registrar en mapas globales:
                mapa_a_primo[orig]   → primo  (para resolver siguiente de pp)
                mapa_a_pp[orig]      → pp     (para resolver anterior de primo)
              Anotar los grupos para la segunda pasada.

    PASADA 2: con los mapas globales completos, asignar siguiente y anterior
              resolviendo cualquier referencia a aristas ya eliminadas.
    """
    # Mapas globales de reemplazo (nombre_orig → objeto_activo)
    mapa_a_primo = {}   # arista eliminada → su primo  (lo que llega  a X)
    mapa_a_pp    = {}   # arista eliminada → su pp     (lo que sale de X)

    todos_los_grupos = []   # lista de (grupos_ordenados, sig_orig[], ant_orig[])

    # ── PASADA 1: crear vértices, primas y pp; registrar mapas ──────────
    ya_procesadas = set()

    for clave, info in intersecciones_reales.items():
        xi, yi = info["punto"].x, info["punto"].y

        semiaristas_orig = []
        vistas = set()
        for seg in info["segmentos"]:
            for arista in mapa_seg_arista.get(id(seg), []):
                for a in (arista, arista.antiarista):
                    if a and a.nombre not in vistas and a.nombre not in ya_procesadas:
                        semiaristas_orig.append(a)
                        vistas.add(a.nombre)

        # Filtrar las que ya tienen origen en X
        semiaristas_orig = [
            a for a in semiaristas_orig
            if not (math.isclose(a.origen.pt.x, xi, abs_tol=1e-7) and
                    math.isclose(a.origen.pt.y, yi, abs_tol=1e-7))
        ]

        if not semiaristas_orig:
            continue

        # Registrar como procesadas
        for a in semiaristas_orig:
            ya_procesadas.add(a.nombre)
    for clave, info in intersecciones_reales.items():
        xi, yi = info["punto"].x, info["punto"].y

        # Reunir semiaristas únicas
        semiaristas_orig = []
        vistas = set()
        for seg in info["segmentos"]:
            for arista in mapa_seg_arista.get(id(seg), []):
                for a in (arista, arista.antiarista):
                    if a and a.nombre not in vistas:
                        semiaristas_orig.append(a)
                        vistas.add(a.nombre)

        if not semiaristas_orig:
            continue

        # Filtrar semi-aristas que ya tienen su origen en X
        semiaristas_orig = [
            a for a in semiaristas_orig
            if not (math.isclose(a.origen.pt.x, xi, abs_tol=1e-7) and
                    math.isclose(a.origen.pt.y, yi, abs_tol=1e-7))
        ]

        if not semiaristas_orig:
            continue
        print(f"\n  Intersección ({xi:.3f}, {yi:.3f}):")
        for a in semiaristas_orig:
            print(
                f"    {a.nombre} | origen={a.origen.pt} | sig={a.siguiente.nombre if a.siguiente else None} | ant={a.anterior.nombre if a.anterior else None}")

        # Crear vértice X
        vx = None
        for v in vertices_g.values():
            if math.isclose(v.pt.x, xi, abs_tol=1e-7) and \
                    math.isclose(v.pt.y, yi, abs_tol=1e-7):
                vx = v
                break

        if vx is None:
            _contador_x[0] += 1
            vx = Vertice(f"X{_contador_x[0]}", Punto(xi, yi))
            vertices_g[vx.nombre] = vx

        # Crear primas y pp, calcular ángulo polar de cada primo
        grupos = []
        for a in semiaristas_orig:
            _contador_a[0] += 1; nom_p  = f"{a.nombre}_p{_contador_a[0]}"
            _contador_a[0] += 1; nom_pp = f"{a.nombre}_pp{_contador_a[0]}"

            a_p  = Arista(nom_p);  a_p.cara  = a.cara
            a_pp = Arista(nom_pp); a_pp.cara = a.cara
            aristas_g[nom_p]  = a_p
            aristas_g[nom_pp] = a_pp

            a_p.origen  = a.origen
            a_pp.origen = vx

            ox, oy = a.origen.pt.x, a.origen.pt.y
            angulo = math.atan2(oy - yi, ox - xi)

            grupos.append({
                "orig":     a,
                "p":        a_p,
                "pp":       a_pp,
                "angulo":   angulo,
                "sig_orig": a.siguiente,   # snapshot antes de eliminar
                "ant_orig": a.anterior,
            })

            # Registrar en mapas globales
            mapa_a_primo[a.nombre] = a_p
            mapa_a_pp   [a.nombre] = a_pp

        # Ordenar CCW y asignar antiaristas (solo dependen del grupo local)
        grupos.sort(key=lambda g: g["angulo"])
        mapa_local = {g["orig"].nombre: g for g in grupos}

        for g in grupos:
            pareja_orig = g["orig"].antiarista
            g_pareja    = mapa_local.get(pareja_orig.nombre) if pareja_orig else None
            if g_pareja:
                g["p"].antiarista  = g_pareja["pp"]
                g["pp"].antiarista = g_pareja["p"]
            else:
                g["p"].antiarista  = g["pp"]
                g["pp"].antiarista = g["p"]

        # Arista incidente de X
        vx.arista_adyacente = grupos[0]["pp"]

        # Actualizar arista_adyacente de vértices de origen
        for g in grupos:
            v_orig = g["p"].origen
            if v_orig and v_orig.arista_adyacente and \
               v_orig.arista_adyacente.nombre == g["orig"].nombre:
                v_orig.arista_adyacente = g["p"]

        # Eliminar originales del dict de aristas
        for g in grupos:
            aristas_g.pop(g["orig"].nombre, None)



        todos_los_grupos.append(grupos)

    # ── PASADA 2: asignar siguiente y anterior con mapas globales ───────
    for grupos in todos_los_grupos:
        N = len(grupos)
        for i, g in enumerate(grupos):
            sig_g = grupos[(i + 1) % N]
            ant_g = grupos[(i - 1) % N]

            # siguiente(primo) = pp del siguiente grupo en CCW [§5.3]
            g["p"].siguiente = sig_g["pp"]

            # siguiente(pp) = siguiente original, resolviendo si fue eliminado [§5.3]
            sig_orig = g["sig_orig"]
            if sig_orig and sig_orig.nombre in mapa_a_primo:
                g["pp"].siguiente = mapa_a_primo[sig_orig.nombre]
            elif sig_orig:
                g["pp"].siguiente = sig_orig
            else:
                g["pp"].siguiente = sig_g["pp"]

            # anterior(primo) = anterior original, resolviendo si fue eliminado [§5.4]
            ant_orig = g["ant_orig"]
            if ant_orig and ant_orig.nombre in mapa_a_pp:
                g["p"].anterior = mapa_a_pp[ant_orig.nombre]
            elif ant_orig:
                g["p"].anterior = ant_orig
            else:
                g["p"].anterior = g["pp"]

            # anterior(pp) = primo del grupo anterior en CCW [§5.4]
            g["pp"].anterior = ant_g["p"]

    # ── PASADA 3: limpieza global ────────────────────────────────────────
    # Aristas que NO fueron subdivididas pueden tener siguiente/anterior
    # apuntando a aristas que SÍ fueron eliminadas en otra intersección.
    # Recorremos todas las aristas supervivientes y resolvemos esas refs.
    for a in list(aristas_g.values()):
        if a.siguiente and a.siguiente.nombre in mapa_a_primo:
            a.siguiente = mapa_a_primo[a.siguiente.nombre]
        if a.anterior and a.anterior.nombre in mapa_a_pp:
            a.anterior = mapa_a_pp[a.anterior.nombre]

    print(f"  {_contador_x[0]} vértices X creados")
    print(f"  DCEL resultante: {len(vertices_g)} vértices, {len(aristas_g)} aristas")


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":

    archivos = ["./layers/layer0" + str(i) for i in range(1, 6)]

    # ── Paso 1 ───────────────────────────────────────────────
    print("── Paso 1: Intersecciones (Bentley-Ottmann) ──")
    segmentos_globales = extraer_todos_los_segmentos(archivos)
    diccionario_intersecciones = encuentraIntersecciones(segmentos_globales)
    print(f"  {len(segmentos_globales)} segmentos, "
          f"{len(diccionario_intersecciones)} puntos de intersección")

    for clave, info in diccionario_intersecciones.items():
        punto = info["punto"]
        for seg in info["segmentos"]:
            if not punto.comparar(seg.p1) and not punto.comparar(seg.p2):
                if punto not in seg.intersecciones:
                    seg.intersecciones.append(punto)

    # ── Paso 2 ───────────────────────────────────────────────
    print("\n── Paso 2: DCEL unificada ──")
    vertices_g, aristas_g, caras_g = cargar_dcel_unificada(archivos)
    mapa_seg_arista       = mapear_segmentos_a_aristas(segmentos_globales, aristas_g)
    intersecciones_reales = filtrar_intersecciones_entre_capas(diccionario_intersecciones)

    # ── Paso 3 ───────────────────────────────────────────────
    print("\n── Paso 3: Subdivisión de aristas (dos pasadas) ──")
    subdividir_todas(intersecciones_reales, mapa_seg_arista, vertices_g, aristas_g)

    # ── Paso 4 ───────────────────────────────────────────────

    def recorrer_ciclo(arista_inicio):
        puntos = []
        actual = arista_inicio
        visitadas = set()
        while actual and actual.nombre not in visitadas:
            visitadas.add(actual.nombre)
            puntos.append((actual.origen.pt.x, actual.origen.pt.y))
            actual = actual.siguiente
        return puntos

    print("\n── Paso 4: Rearme de caras ──")
    constructor = ConstructorCaras(vertices_g, aristas_g)
    caras_g = constructor.construir()
    print("\n  Diagnóstico de caras:")
    for nombre, cara in caras_g.items():
        if nombre == "f_infinita":
            continue
        if cara.aristas_externas is None:
            print(f"  ⚠️  {nombre}: aristas_externas = None")
            continue
        puntos = recorrer_ciclo(cara.aristas_externas)
        xs = [p[0] for p in puntos]
        ys = [p[1] for p in puntos]
        print(f"  {nombre}: {len(puntos)} vértices | "
              f"x=[{min(xs):.1f},{max(xs):.1f}] "
              f"y=[{min(ys):.1f},{max(ys):.1f}]")



    print("\n  Detalle de f0:")
    cara = caras_g["f0"]
    actual = cara.aristas_externas
    visitadas = set()
    while actual and actual.nombre not in visitadas:
        visitadas.add(actual.nombre)
        print(f"    {actual.nombre} | origen={actual.origen.pt} | sig={actual.siguiente.nombre}")
        actual = actual.siguiente

    print(f"  {len(caras_g)} caras encontradas (incluyendo f_infinita)")

    # Verificación de integridad
    print("\n  Verificación de integridad:")
    errores = 0
    for nombre, a in aristas_g.items():
        for campo, val in [("origen",     a.origen),
                           ("antiarista", a.antiarista),
                           ("siguiente",  a.siguiente),
                           ("anterior",   a.anterior)]:
            if val is None:
                print(f"    ⚠️  {nombre}.{campo} = None"); errores += 1
            elif campo != "origen" and hasattr(val, 'nombre') and val.nombre not in aristas_g:
                print(f"    ⚠️  {nombre}.{campo} → '{val.nombre}' no existe"); errores += 1
    if errores == 0:
        print("    ✅ Todas las referencias son válidas")
    else:
        print(f"    ❌ {errores} referencias inválidas")

    # ── Visualización ────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 10))
    colores = ['blue', 'green', 'orange', 'purple', 'brown']
    colores_caras = cm.Set3.colors


    # Relleno de caras (zorder=1 para que quede debajo de todo)



    for i, (nombre, cara) in enumerate(caras_g.items()):
        if nombre == "f_infinita" or cara.aristas_externas is None:
            continue
        puntos = recorrer_ciclo(cara.aristas_externas)
        if len(puntos) < 3:
            continue
        xs = [p[0] for p in puntos]
        ys = [p[1] for p in puntos]
        color = colores_caras[i % len(colores_caras)] + (0.3,)
        ax.fill(xs, ys, color=color, zorder=1)
        cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
        ax.text(cx, cy, nombre, fontsize=7, ha='center', va='center', zorder=2)

    # Segmentos originales por capa
    for seg in segmentos_globales:
        ax.plot([seg.p1.x, seg.p2.x], [seg.p1.y, seg.p2.y],
                color=colores[seg.layer_id % len(colores)],
                linewidth=2.5, alpha=0.18, label=f"Layer {seg.layer_id + 1}")

    # Aristas de la DCEL
    dibujadas = set()
    for nombre, a in aristas_g.items():
        if nombre in dibujadas: continue
        if not a.antiarista or a.antiarista.nombre not in aristas_g: continue
        p1 = a.origen.pt
        p2 = a.antiarista.origen.pt
        ax.plot([p1.x, p2.x], [p1.y, p2.y], 'k-', linewidth=1.2, alpha=0.9, zorder=3)
        dibujadas.add(nombre)
        dibujadas.add(a.antiarista.nombre)

    # Vértices
    for v in vertices_g.values():
        es_x = v.nombre.startswith('X')
        ax.plot(v.pt.x, v.pt.y, 'o',
                color='red' if es_x else 'black',
                markersize=6 if es_x else 3, zorder=5)
        if es_x:
            ax.text(v.pt.x + 0.1, v.pt.y + 0.1, v.nombre, fontsize=6, color='red')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels, handles)).values(),
              dict(zip(labels, handles)).keys(), loc='upper right')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.4)
    plt.title("Paso 4 – DCEL con caras")
    plt.tight_layout()
    plt.show()
