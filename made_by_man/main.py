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
Algoritmo fusión de listas de aristas / fusión de caras — MAPOVERLAY (De Berg §2.3)

ESTADO ACTUAL → Paso 4 completo:
  ✅ Paso 1 – Intersecciones (Bentley-Ottmann)
  ✅ Paso 2 – DCEL unificada con prefijos por capa
  ✅ Paso 3 – Subdivisión con pre-fragmentación por arista:
              Cada arista que pasa por N intersecciones se divide en N+1
              fragmentos de una sola vez (ordenados por t), eliminando el
              bug donde mapa_a_primo solo guardaba el último primo.
              Las antiaristas se asignan emparejando fragmentos en sentido
              inverso. La Fase C usa frag_de para resolver correctamente
              las referencias de aristas no subdivididas.
  ✅ Paso 4 – Rearme de caras (De Berg §2.3 pasos 4-7):
              4a. Extraer ciclos recorriendo siguiente
              4b. Clasificar exterior/interior (producto cruzado en vértice más izquierdo)
              4c. Construir grafo G: hole → ciclo a su izquierda
              4d. Componentes conexas de G → una cara por componente
              4e. Asignar campo .cara a cada semiarista
"""

from Figuras import *
from Graficacion import *
from ArbolBarrido import *
from Segmentacion import *
import heapq
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    Prefixa cada nombre con 'L{i}_' para evitar colisiones entre capas
    (p.ej. 'p1' en layer01 y 'p1' en layer03 son vértices distintos).
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
# PASO 3 – Subdivisión con pre-fragmentación por arista
# ─────────────────────────────────────────────────────────────

_contador_x = [0]
_contador_a = [0]


def _pts_eq(p, x, y):
    return math.isclose(p.x, x, abs_tol=1e-7) and math.isclose(p.y, y, abs_tol=1e-7)


def _param_t(orig_pt, dest_pt, xi, yi):
    """Parámetro t ∈ [0,1] del punto (xi,yi) sobre el segmento orig→dest."""
    dx = dest_pt.x - orig_pt.x
    dy = dest_pt.y - orig_pt.y
    if abs(dx) > abs(dy):
        return (xi - orig_pt.x) / dx if abs(dx) > 1e-12 else 0.0
    return (yi - orig_pt.y) / dy if abs(dy) > 1e-12 else 0.0


def subdividir_todas(intersecciones_reales, mapa_seg_arista, vertices_g, aristas_g):
    """
    Subdivide todas las intersecciones en tres fases.

    FASE A – Pre-fragmentación por arista:
      Para cada arista que pasa por N intersecciones, crea N+1 fragmentos
      de una vez, ordenados por parámetro t. Esto evita el bug anterior donde
      una arista con múltiples intersecciones producía primos incorrectos.
      Los pares (a, anti_a) se procesan juntos para asignar antiaristas
      correctamente: frag_a[i].antiarista = frag_anti[N-i].

    FASE B – Grupos polares por intersección:
      Para cada Xi, reúne los pares (primo, pp) de cada semiarista que pasa
      por ella, aplica ordenamiento polar CCW y asigna siguiente/anterior/
      antiarista según las reglas §5 del pseudocódigo.

    FASE C – Limpieza global:
      Resuelve referencias de aristas supervivientes (no subdivididas) que
      aún apuntan a aristas eliminadas, usando frag_de para garantizar que
      se resuelve al fragmento correcto (frag[0] para .siguiente, frag[-1]
      para .anterior).
    """
    # ── Paso previo: crear todos los vértices X ───────────────────────────
    vertices_x = {}
    for clave, info in intersecciones_reales.items():
        xi, yi = info["punto"].x, info["punto"].y
        _contador_x[0] += 1
        vx = Vertice(f"X{_contador_x[0]}", Punto(xi, yi))
        vertices_g[vx.nombre] = vx
        vertices_x[clave] = vx

    # ── FASE A ────────────────────────────────────────────────────────────
    # Recopilar intersecciones interiores por cada semiarista.
    # Usamos un set de claves (nombre, clave) para evitar duplicados
    # (el mismo punto puede aparecer múltiples veces si varios segmentos
    # geométricos mapean a la misma semiarista DCEL).
    arista_a_puntos = {}   # nombre → [(t, clave, vx), ...]
    arista_puntos_vistos = set()   # (nombre, clave) ya procesado
    for clave, info in intersecciones_reales.items():
        xi, yi = info["punto"].x, info["punto"].y
        for seg in info["segmentos"]:
            for a in mapa_seg_arista.get(id(seg), []):
                for ha in (a, a.antiarista):
                    if ha is None: continue
                    dedup_key = (ha.nombre, clave)
                    if dedup_key in arista_puntos_vistos:
                        continue
                    arista_puntos_vistos.add(dedup_key)
                    orig = ha.origen.pt
                    dest = ha.antiarista.origen.pt
                    if _pts_eq(orig, xi, yi) or _pts_eq(dest, xi, yi):
                        continue   # punto en extremo, no interior
                    t = _param_t(orig, dest, xi, yi)
                    arista_a_puntos.setdefault(ha.nombre, []).append(
                        (t, clave, vertices_x[clave]))

    primo_en = {}   # (nombre_arista, clave) → fragmento primo
    pp_en    = {}   # (nombre_arista, clave) → fragmento pp
    frag_de  = {}   # nombre_arista → [f0, f1, ..., fN]
    aristas_eliminadas = set()
    pares_procesados   = set()

    def _crear_frags(orig_arista, puntos_lista):
        """Crea N+1 fragmentos para orig_arista con cortes en puntos_lista (ya ordenados por t)."""
        N = len(puntos_lista)
        frags = []
        for _ in range(N + 1):
            _contador_a[0] += 1
            f = Arista(f"{orig_arista.nombre}_f{_contador_a[0]}")
            f.cara = orig_arista.cara
            aristas_g[f.nombre] = f
            frags.append(f)
        frags[0].origen = orig_arista.origen
        for j in range(1, N + 1):
            frags[j].origen = puntos_lista[j - 1][2]
        # Extremos heredan del original
        frags[0].anterior    = orig_arista.anterior
        frags[N].siguiente   = orig_arista.siguiente
        # Provisionales internos (sobreescritos en Fase B)
        for j in range(N + 1):
            if j < N: frags[j].siguiente = frags[j + 1]
            if j > 0: frags[j].anterior  = frags[j - 1]
        return frags

    for nombre_a in list(arista_a_puntos.keys()):
        if nombre_a in pares_procesados or nombre_a not in aristas_g:
            continue
        a      = aristas_g[nombre_a]
        anti_a = a.antiarista
        nombre_anti = anti_a.nombre if anti_a else None

        puntos_a    = sorted(arista_a_puntos.get(nombre_a,    []), key=lambda x: x[0])
        puntos_anti = sorted(arista_a_puntos.get(nombre_anti, []), key=lambda x: x[0]) \
                      if nombre_anti else []

        frags_a    = _crear_frags(a,    puntos_a)
        frags_anti = _crear_frags(anti_a, puntos_anti) \
                     if anti_a and puntos_anti else []

        # Antiaristas: frag_a[i] ↔ frag_anti[N-i]  (sentidos opuestos)
        N_a = len(puntos_a); N_anti = len(puntos_anti)
        if frags_anti and N_a == N_anti:
            for j in range(N_a + 1):
                frags_a[j].antiarista              = frags_anti[N_anti - j]
                frags_anti[N_anti - j].antiarista  = frags_a[j]
        else:
            for f in frags_a:    f.antiarista = anti_a
            for f in frags_anti: f.antiarista = a

        for j, (t, clave, vx) in enumerate(puntos_a):
            primo_en[(nombre_a,    clave)] = frags_a[j]
            pp_en   [(nombre_a,    clave)] = frags_a[j + 1]
        for j, (t, clave, vx) in enumerate(puntos_anti):
            primo_en[(nombre_anti, clave)] = frags_anti[j]
            pp_en   [(nombre_anti, clave)] = frags_anti[j + 1]

        frag_de[nombre_a] = frags_a
        if frags_anti: frag_de[nombre_anti] = frags_anti

        aristas_eliminadas.add(nombre_a)
        if nombre_anti and nombre_anti in arista_a_puntos:
            aristas_eliminadas.add(nombre_anti)
            pares_procesados.add(nombre_anti)
        pares_procesados.add(nombre_a)

    for nombre_a in aristas_eliminadas:
        aristas_g.pop(nombre_a, None)

    # ── FASE B ────────────────────────────────────────────────────────────
    for clave, info in intersecciones_reales.items():
        xi, yi = info["punto"].x, info["punto"].y
        vx = vertices_x[clave]

        # Reunir pares (primo, pp) para cada semiarista que pasa por Xi
        pares = []
        vistos_pares = set()
        for seg in info["segmentos"]:
            for a in mapa_seg_arista.get(id(seg), []):
                for ha_nombre in (a.nombre,
                                  a.antiarista.nombre if a.antiarista else None):
                    if ha_nombre is None or ha_nombre in vistos_pares:
                        continue
                    key = (ha_nombre, clave)
                    if key in primo_en:
                        pares.append({
                            "primo":       primo_en[key],
                            "pp":          pp_en   [key],
                            "orig_nombre": ha_nombre,
                        })
                        vistos_pares.add(ha_nombre)

        if not pares:
            continue   # intersección en vértice compartido, sin cortes interiores

        # Ángulo del primo desde X → hacia su origen
        for par in pares:
            ox = par["primo"].origen.pt.x
            oy = par["primo"].origen.pt.y
            par["angulo"] = math.atan2(oy - yi, ox - xi)
        pares.sort(key=lambda p: p["angulo"])
        N = len(pares)

        # Antiaristas entre pares: anti(primo_a) = pp_pareja, anti(pp_a) = primo_pareja
        nombre_a_par = {p["orig_nombre"]: p for p in pares}
        for par in pares:
            pareja_nombre = None
            for seg in info["segmentos"]:
                for a in mapa_seg_arista.get(id(seg), []):
                    if a.nombre == par["orig_nombre"] and a.antiarista:
                        pareja_nombre = a.antiarista.nombre; break
                    if a.antiarista and a.antiarista.nombre == par["orig_nombre"]:
                        pareja_nombre = a.nombre; break
                if pareja_nombre: break
            par_pareja = nombre_a_par.get(pareja_nombre) if pareja_nombre else None
            if par_pareja:
                par["primo"].antiarista = par_pareja["pp"]
                par["pp"].antiarista    = par_pareja["primo"]
            else:
                par["primo"].antiarista = par["pp"]
                par["pp"].antiarista    = par["primo"]

        # Siguiente y anterior (reglas §5.3 y §5.4)
        for i, par in enumerate(pares):
            sig_par = pares[(i + 1) % N]
            ant_par = pares[(i - 1) % N]
            # §5.3  siguiente(primo) = pp del siguiente en CCW
            par["primo"].siguiente = sig_par["pp"]
            # §5.4  anterior(pp)    = primo del anterior en CCW
            par["pp"].anterior     = ant_par["primo"]
            # §5.3  siguiente(pp)   = siguiente(primo) del mismo fragmento
            #       (el pp y el primo son el mismo objeto frag[i]; cuando
            #       la Fase B procesa Xi+1 sobreescribe el sig del primo,
            #       pero aquí en Xi el pp es ese mismo objeto y su sig debe
            #       ser la arista que sigue en el ciclo PASANDO por Xi+1.
            #       Esa arista es precisamente el pp del siguiente en CCW en Xi+1,
            #       que es par["primo"].siguiente (lo que acaba de asignarse).
            #       Pero par["primo"] Y par["pp"] son el mismo objeto frag[i],
            #       así que esto es redundante — el sig ya está asignado arriba.
            # §5.4  anterior(primo) viene de Fase A (hereda del original si es frag[0],
            #       o viene del orden polar si es intermedio — se deja sin tocar aquí
            #       para que el provisional de Fase A lo complete).

        # Arista incidente del vértice X
        vx.arista_adyacente = pares[0]["pp"]

        # Actualizar arista_adyacente de vértices de origen cuyo incidente fue eliminado
        for par in pares:
            v_orig = par["primo"].origen
            if v_orig and v_orig.arista_adyacente and \
               v_orig.arista_adyacente.nombre in aristas_eliminadas:
                v_orig.arista_adyacente = par["primo"]

    # ── FASE C: limpieza global ───────────────────────────────────────────
    # frag_de garantiza el fragmento correcto: frag[0] mismo origen, frag[-1] mismo destino
    mapa_sig = {nombre: frags[0]  for nombre, frags in frag_de.items()}
    mapa_ant = {nombre: frags[-1] for nombre, frags in frag_de.items()}

    for a in list(aristas_g.values()):
        if a.siguiente and a.siguiente.nombre in aristas_eliminadas:
            a.siguiente = mapa_sig.get(a.siguiente.nombre, a.siguiente)
        if a.anterior and a.anterior.nombre in aristas_eliminadas:
            a.anterior  = mapa_ant.get(a.anterior.nombre,  a.anterior)

    print(f"  {_contador_x[0]} vértices X creados")
    print(f"  DCEL resultante: {len(vertices_g)} vértices, {len(aristas_g)} aristas")


# ─────────────────────────────────────────────────────────────
# PASO 4 – Rearme de caras  (De Berg §2.3, pasos 4-7)
# ─────────────────────────────────────────────────────────────

def extraer_ciclos(aristas_g):
    """
    Recorre la DCEL siguiendo .siguiente y extrae todos los ciclos.
    Cada ciclo es una lista de objetos Arista en orden.
    Garantía: cada arista aparece en exactamente un ciclo.
    """
    visitadas = set()
    ciclos    = []
    for nombre_inicio in list(aristas_g.keys()):
        if nombre_inicio in visitadas:
            continue
        ciclo  = []
        actual = aristas_g.get(nombre_inicio)
        while actual and actual.nombre not in visitadas:
            visitadas.add(actual.nombre)
            ciclo.append(actual)
            actual = actual.siguiente
        if ciclo:
            ciclos.append(ciclo)
    return ciclos


def _area_signed(ciclo):
    """Área con signo (Shoelace). Positiva=CCW, negativa=CW."""
    pts = [(a.origen.pt.x, a.origen.pt.y) for a in ciclo]
    n   = len(pts)
    return sum(
        pts[j][0] * pts[(j+1)%n][1] - pts[(j+1)%n][0] * pts[j][1]
        for j in range(n)
    ) / 2.0


def _clasificar_ciclo(ciclo):
    """
    True → ciclo exterior (outer boundary, recorrido CCW, área > 0).
    False → ciclo interior (hole o cara infinita, recorrido CW, área < 0).

    Método primario: área con signo (Shoelace), robusto para todos los casos.
    Fallback para ciclos degenerados (área ≈ 0): producto cruzado en el
    vértice más izquierdo (método De Berg §2.3).
    """
    area = _area_signed(ciclo)
    if abs(area) > 1e-9:
        return area > 0

    # Fallback: producto cruzado en vértice más izquierdo
    idx_min = 0
    for i, a in enumerate(ciclo):
        v  = a.origen.pt
        v0 = ciclo[idx_min].origen.pt
        if v.x < v0.x or (math.isclose(v.x, v0.x, abs_tol=1e-9) and v.y < v0.y):
            idx_min = i
    n        = len(ciclo)
    a_actual = ciclo[idx_min]
    a_ent    = ciclo[(idx_min - 1) % n]
    v_izq    = a_actual.origen.pt
    v_ent    = a_ent.origen.pt
    if a_actual.antiarista is None:
        return False
    v_sal  = a_actual.antiarista.origen.pt
    dx_in  = v_izq.x - v_ent.x;  dy_in  = v_izq.y - v_ent.y
    dx_out = v_sal.x - v_izq.x;  dy_out = v_sal.y - v_izq.y
    return (dx_in * dy_out - dy_in * dx_out) > 0


def _arista_a_la_izquierda(ciclo, aristas_g):
    """
    Para el vértice más izquierdo de un ciclo (que es hole), busca la
    semiarista de la DCEL que está inmediatamente a la izquierda de ese
    vértice a la misma altura y.

    De Berg: 'If e is the half-edge immediately to the left of v, then
    we add an arc between the cycle containing e and the hole cycle.'

    La semiarista buscada es la que:
      - Su segmento geométrico cruza la horizontal y = v_izq.y
      - Su x evaluado en esa y es el máximo estrictamente menor que v_izq.x
      - Está orientada de modo que la cara queda arriba (apunta hacia la derecha),
        es decir, el origen tiene y >= v_izq.y (viene de arriba o del mismo nivel).
    """
    idx_min = 0
    for i, a in enumerate(ciclo):
        v  = a.origen.pt
        v0 = ciclo[idx_min].origen.pt
        if v.x < v0.x or (math.isclose(v.x, v0.x, abs_tol=1e-9) and v.y < v0.y):
            idx_min = i

    v_izq = ciclo[idx_min].origen.pt
    y_ref = v_izq.y

    mejor   = None
    mejor_x = -math.inf

    for a in aristas_g.values():
        if a.antiarista is None: continue
        p1 = a.origen.pt
        p2 = a.antiarista.origen.pt

        y_min = min(p1.y, p2.y)
        y_max = max(p1.y, p2.y)

        # El segmento debe cruzar la línea horizontal y = y_ref
        if y_ref < y_min - 1e-9 or y_ref > y_max + 1e-9:
            continue

        # x del segmento evaluado en y_ref
        if abs(p2.y - p1.y) < 1e-9:
            x_seg = min(p1.x, p2.x)
        else:
            t = (y_ref - p1.y) / (p2.y - p1.y)
            t = max(0.0, min(1.0, t))
            x_seg = p1.x + t * (p2.x - p1.x)

        if x_seg >= v_izq.x - 1e-9:
            continue   # está a la derecha o en el mismo punto

        # Solo la semiarista que apunta hacia la derecha (la cara queda arriba):
        # su origen tiene y >= y_ref (viene de arriba o del mismo nivel)
        if p1.y < y_ref - 1e-9:
            continue   # origen está debajo de y_ref → apunta hacia arriba, cara abajo

        if x_seg > mejor_x:
            mejor_x = x_seg
            mejor   = a

    return mejor


def rearmar_caras(ciclos, aristas_g):
    """
    Implementa De Berg §2.3 pasos 4-7:

    1. Clasifica cada ciclo como exterior o interior.
    2. Para cada ciclo interior (hole), busca la semiarista inmediatamente
       a la izquierda de su vértice más izquierdo → arco en el grafo G.
       Si no hay semiarista a la izquierda → arco a la cara no acotada (∞).
    3. Cada componente conexa del grafo G = una cara del overlay.
       - La componente tiene exactamente un ciclo exterior (su outer boundary).
       - El resto de ciclos de la componente son holes (inner boundaries).
    4. Crea objetos Cara y asigna .cara a todas las semiaristas.

    Retorna el diccionario de caras nuevas.
    """
    # ── 4a. Clasificar ciclos ─────────────────────────────────────────────
    exteriores = [c for c in ciclos if     _clasificar_ciclo(c)]
    interiores  = [c for c in ciclos if not _clasificar_ciclo(c)]
    print(f"  Ciclos exteriores: {len(exteriores)}, holes: {len(interiores)}")

    # Índice: nombre_de_cualquier_arista_del_ciclo → índice_del_ciclo
    id_de_ciclo = {}
    for i, ciclo in enumerate(ciclos):
        for a in ciclo:
            id_de_ciclo[a.nombre] = i

    N_ciclos = len(ciclos)
    INF      = N_ciclos   # nodo especial: cara no acotada

    # ── 4b. Construir grafo G ─────────────────────────────────────────────
    # grafo[i] = j significa que el ciclo i (hole) se conecta al ciclo j
    grafo = {}
    for ciclo in interiores:
        idx_hole = id_de_ciclo[ciclo[0].nombre]
        vecina   = _arista_a_la_izquierda(ciclo, aristas_g)
        if vecina and vecina.nombre in id_de_ciclo:
            grafo[idx_hole] = id_de_ciclo[vecina.nombre]
        else:
            grafo[idx_hole] = INF   # hole sin nada a la izquierda → cara infinita

    # ── 4c. Componentes conexas de G (BFS/union-find) ────────────────────
    # Representación: padre[i] = índice del padre en el grafo
    # Encontrar la raíz de cada ciclo transitivamente
    raiz_cache = {}

    def raiz(i):
        if i not in raiz_cache:
            if i == INF or i not in grafo:
                raiz_cache[i] = i
            else:
                raiz_cache[i] = raiz(grafo[i])
        return raiz_cache[i]

    # Agrupar ciclos por raíz
    componentes = {}   # raiz_idx → {"ext": ciclo, "holes": [ciclo, ...]}
    for i, ciclo in enumerate(ciclos):
        r = raiz(i)
        if r not in componentes:
            componentes[r] = {"ext": None, "holes": []}
        if ciclo in exteriores:
            componentes[r]["ext"] = ciclo
        else:
            componentes[r]["holes"].append(ciclo)

    if INF not in componentes:
        componentes[INF] = {"ext": None, "holes": []}

    # ── 4d. Crear objetos Cara ────────────────────────────────────────────
    nuevas_caras = {}
    cara_idx     = [0]

    for r, comp in componentes.items():
        cara_idx[0] += 1
        nombre_cara = "C_inf" if r == INF else f"C{cara_idx[0]}"
        cara = Cara(nombre_cara)

        ciclo_ext  = comp["ext"]
        ciclos_int = comp["holes"]

        # aristas_externas: semiarista del outer boundary (borde de la cara)
        cara.aristas_externas = ciclo_ext[0] if ciclo_ext else None

        # aristas_internas: una semiarista representativa de cada hole
        cara.aristas_internas = [c[0] for c in ciclos_int]

        nuevas_caras[nombre_cara] = cara

        # ── 4e. Asignar .cara a cada semiarista ───────────────────────────
        todos_ciclos = ([ciclo_ext] if ciclo_ext else []) + ciclos_int
        for ciclo in todos_ciclos:
            for a in ciclo:
                a.cara = cara

    return nuevas_caras


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
    print("\n── Paso 3: Subdivisión de aristas ──")
    subdividir_todas(intersecciones_reales, mapa_seg_arista, vertices_g, aristas_g)

    # Verificación de integridad
    print("\n  Verificación de integridad DCEL:")
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

    # Verificar que todos los ciclos son cerrados (sin estructuras ρ)
    rho_count = 0
    for nombre in aristas_g:
        a = aristas_g[nombre]
        local = {}
        act = a
        for _ in range(len(aristas_g) + 2):
            if act is None or act.nombre in local:
                if act and act.nombre != a.nombre:
                    rho_count += 1
                break
            local[act.nombre] = True
            act = act.siguiente

    if errores == 0 and rho_count == 0:
        print("    ✅ DCEL íntegra, todos los ciclos cerrados")
    else:
        if errores:   print(f"    ❌ {errores} referencias inválidas")
        if rho_count: print(f"    ❌ {rho_count} cadenas con estructura ρ (ciclos no cerrados)")

    # ── Paso 4 ───────────────────────────────────────────────
    print("\n── Paso 4: Rearme de caras (De Berg §2.3) ──")
    ciclos = extraer_ciclos(aristas_g)
    print(f"  Ciclos extraídos: {len(ciclos)}  "
          f"(tamaños: {sorted(len(c) for c in ciclos)})")

    caras_g = rearmar_caras(ciclos, aristas_g)
    print(f"\n  Caras construidas: {len(caras_g)}")
    for nombre, cara in sorted(caras_g.items()):
        ext  = cara.aristas_externas.nombre if cara.aristas_externas else "—"
        ints = [a.nombre for a in cara.aristas_internas]
        print(f"    {nombre}: exterior={ext}  holes={ints}")

    # ── Visualización ─────────────────────────────────────────
    import random; random.seed(42)
    colores_cara = {}
    for nombre in caras_g:
        if nombre != "C_inf":
            colores_cara[nombre] = (random.random(), random.random(), random.random())

    fig, ax = plt.subplots(figsize=(12, 12))

    # Rellenar cada cara acotada
    for nombre, cara in caras_g.items():
        if nombre == "C_inf" or cara.aristas_externas is None:
            continue
        color = colores_cara.get(nombre, (0.8, 0.8, 0.8))
        pts = []
        actual = cara.aristas_externas
        for _ in range(500):
            pts.append((actual.origen.pt.x, actual.origen.pt.y))
            actual = actual.siguiente
            if actual is None or actual is cara.aristas_externas:
                break
        if len(pts) >= 3:
            xs, ys = zip(*pts)
            ax.fill(xs, ys, color=color, alpha=0.25, zorder=1)

    # Aristas
    dibujadas = set()
    for nombre, a in aristas_g.items():
        if nombre in dibujadas: continue
        if not a.antiarista or a.antiarista.nombre not in aristas_g: continue
        p1 = a.origen.pt; p2 = a.antiarista.origen.pt
        ax.plot([p1.x, p2.x], [p1.y, p2.y], 'k-', linewidth=1.0, alpha=0.85, zorder=3)
        dibujadas.add(nombre); dibujadas.add(a.antiarista.nombre)

    # Vértices
    for v in vertices_g.values():
        es_x = v.nombre.startswith('X')
        ax.plot(v.pt.x, v.pt.y, 'o',
                color='red' if es_x else 'black',
                markersize=5 if es_x else 3, zorder=5)
        if es_x:
            ax.text(v.pt.x + 0.08, v.pt.y + 0.08, v.nombre, fontsize=5, color='red')

    # Leyenda de caras
    patches = [mpatches.Patch(color=colores_cara[n], label=n)
               for n in sorted(colores_cara)]
    ax.legend(handles=patches, loc='upper right', fontsize=7,
              title="Caras", title_fontsize=8)

    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.35)
    plt.title(f"Paso 4 – {len(caras_g)} caras reconstruidas  "
              f"({len(ciclos)} ciclos, {len(aristas_g)} semiaristas)")
    plt.tight_layout()
    plt.savefig("paso4_caras.png", dpi=150)
    print("\n  Gráfica guardada en paso4_caras.png")
    plt.show()
