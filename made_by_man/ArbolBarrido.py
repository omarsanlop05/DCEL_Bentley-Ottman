# -*- coding: utf-8 -*-
"""
ArbolBarrido.py  –  versión iterativa + AVL balanceado
Árbol AVL para la línea de barrido del algoritmo Bentley-Ottmann.

La línea de barrido es HORIZONTAL y avanza de arriba hacia abajo (Y descendente).
El orden de los segmentos se determina por su coordenada X evaluada en la posición
actual de la línea de barrido (sweep_y).

Convención de segmentos:
    p1 = extremo superior (mayor Y)
    p2 = extremo inferior (menor Y)

Esta versión usa un AVL con operaciones 100% iterativas para soportar
conjuntos de entrada muy grandes (>100k segmentos) sin riesgo de
RecursionError.
"""


# ──────────────────────────────────────────────────────────────────────────────
# Nodo AVL
# ──────────────────────────────────────────────────────────────────────────────

class _Nodo:
    __slots__ = ("segmento", "izq", "der", "padre", "altura")

    def __init__(self, segmento):
        self.segmento = segmento
        self.izq   = None
        self.der   = None
        self.padre = None
        self.altura = 1


def _altura(n):
    return n.altura if n else 0


def _actualizar_altura(n):
    if n:
        n.altura = 1 + max(_altura(n.izq), _altura(n.der))


def _factor(n):
    return _altura(n.izq) - _altura(n.der) if n else 0


# ──────────────────────────────────────────────────────────────────────────────
# Árbol AVL iterativo
# ──────────────────────────────────────────────────────────────────────────────

class ArbolBarrido:
    """
    AVL que ordena segmentos por su coordenada X en sweep_y.
    Todas las operaciones (insertar, eliminar, buscar, en_orden, vecinos)
    son iterativas: no hay riesgo de RecursionError con millones de nodos.

    Atributos públicos
    ------------------
    sweep_y : float
        Posición actual de la línea de barrido horizontal.
        Debe actualizarse antes de llamar a insertar / eliminar / vecinos.
    """

    def __init__(self):
        self._raiz  = None
        self.sweep_y = 0.0

    # ── utilidades de clave ───────────────────────────────────────────────────

        # En ArbolBarrido.py, modifica estos dos métodos:

    def _x_en_sweep(self, seg):
        y1, x1 = seg.p1.y, seg.p1.x
        y2, x2 = seg.p2.y, seg.p2.x

        # Si el segmento es horizontal, devolvemos el x menor para estabilidad
        if abs(y2 - y1) < 1e-9:
            return min(x1, x2)

        # Interpolación con límite de seguridad
        t = (self.sweep_y - y1) / (y2 - y1)
        # Forzamos que t esté en el rango [0, 1] por precisión numérica
        t = max(0.0, min(1.0, t))
        return x1 + t * (x2 - x1)

    def _clave(self, seg):
        x = self._x_en_sweep(seg)
        dy = seg.p2.y - seg.p1.y

        # Mejoramos el cálculo de la pendiente para evitar inf
        if abs(dy) < 1e-9:
            slope = 1e12  # Un valor muy alto pero finito
        else:
            slope = (seg.p2.x - seg.p1.x) / dy

        # Redondeamos ligeramente la X para que el árbol sea tolerante a errores de float
        return (round(x, 9), round(slope, 9))

    def _mismo_seg(self, a, b):
        return (a.p1 == b.p1 and a.p2 == b.p2) or \
               (a.p1 == b.p2 and a.p2 == b.p1)

    # ── rotaciones AVL ────────────────────────────────────────────────────────

    def _rotar_der(self, y):
        x   = y.izq
        T2  = x.der

        x.der  = y
        y.izq  = T2

        # actualizar padres
        x.padre = y.padre
        y.padre = x
        if T2:
            T2.padre = y

        if x.padre is None:
            self._raiz = x
        elif x.padre.izq is y:
            x.padre.izq = x
        else:
            x.padre.der = x

        _actualizar_altura(y)
        _actualizar_altura(x)
        return x

    def _rotar_izq(self, x):
        y   = x.der
        T2  = y.izq

        y.izq  = x
        x.der  = T2

        y.padre = x.padre
        x.padre = y
        if T2:
            T2.padre = x

        if y.padre is None:
            self._raiz = y
        elif y.padre.izq is x:
            y.padre.izq = y
        else:
            y.padre.der = y

        _actualizar_altura(x)
        _actualizar_altura(y)
        return y

    def _balancear(self, nodo):
        """Sube desde nodo hacia la raíz rebalanceando."""
        n = nodo
        while n:
            _actualizar_altura(n)
            fb = _factor(n)

            if fb > 1:                          # pesado a la izquierda
                if _factor(n.izq) < 0:
                    self._rotar_izq(n.izq)
                n = self._rotar_der(n)
            elif fb < -1:                       # pesado a la derecha
                if _factor(n.der) > 0:
                    self._rotar_der(n.der)
                n = self._rotar_izq(n)

            n = n.padre

    # ── insertar (iterativo) ──────────────────────────────────────────────────

    def insertar(self, segmento):
        nuevo = _Nodo(segmento)
        if self._raiz is None:
            self._raiz = nuevo
            return

        actual = self._raiz
        while True:
            if self._mismo_seg(segmento, actual.segmento):
                return                          # ya existe
            if self._clave(segmento) < self._clave(actual.segmento):
                if actual.izq is None:
                    actual.izq = nuevo
                    nuevo.padre = actual
                    break
                actual = actual.izq
            else:
                if actual.der is None:
                    actual.der = nuevo
                    nuevo.padre = actual
                    break
                actual = actual.der

        self._balancear(nuevo.padre)

    # ── eliminar (iterativo) ──────────────────────────────────────────────────

    def eliminar(self, segmento):
        nodo = self._buscar_nodo(segmento)
        if nodo is None:
            return
        self._eliminar_nodo(nodo)

    def _eliminar_nodo(self, nodo):
        # Caso 1: hoja
        if nodo.izq is None and nodo.der is None:
            padre = nodo.padre
            self._reemplazar_en_padre(nodo, None)
            self._balancear(padre)

        # Caso 2: solo hijo derecho
        elif nodo.izq is None:
            padre = nodo.padre
            self._reemplazar_en_padre(nodo, nodo.der)
            nodo.der.padre = padre
            self._balancear(padre)

        # Caso 3: solo hijo izquierdo
        elif nodo.der is None:
            padre = nodo.padre
            self._reemplazar_en_padre(nodo, nodo.izq)
            nodo.izq.padre = padre
            self._balancear(padre)

        # Caso 4: dos hijos → sucesor in-order
        else:
            sucesor = nodo.der
            while sucesor.izq:
                sucesor = sucesor.izq
            nodo.segmento = sucesor.segmento   # copiar dato
            self._eliminar_nodo(sucesor)        # eliminar sucesor (≤1 hijo)

    def _reemplazar_en_padre(self, nodo, reemplazo):
        if nodo.padre is None:
            self._raiz = reemplazo
        elif nodo.padre.izq is nodo:
            nodo.padre.izq = reemplazo
        else:
            nodo.padre.der = reemplazo
        if reemplazo is not None:
            reemplazo.padre = nodo.padre

    # ── búsqueda (iterativa) ──────────────────────────────────────────────────

    def buscar(self, segmento):
        return self._buscar_nodo(segmento) is not None

    def _buscar_nodo(self, segmento):
        actual = self._raiz
        while actual:
            if self._mismo_seg(segmento, actual.segmento):
                return actual
            if self._clave(segmento) < self._clave(actual.segmento):
                actual = actual.izq
            else:
                actual = actual.der
        return None

    # ── vecinos (iterativo) ───────────────────────────────────────────────────

    def vecinos(self, segmento):
        """Devuelve (izquierdo, derecho) del segmento en el orden actual."""
        izq = None
        der = None
        nodo = self._raiz

        while nodo:
            if self._mismo_seg(segmento, nodo.segmento):
                # predecesor: máximo del sub-árbol izquierdo
                if nodo.izq:
                    tmp = nodo.izq
                    while tmp.der:
                        tmp = tmp.der
                    izq = tmp.segmento
                # sucesor: mínimo del sub-árbol derecho
                if nodo.der:
                    tmp = nodo.der
                    while tmp.izq:
                        tmp = tmp.izq
                    der = tmp.segmento
                break
            elif self._clave(segmento) < self._clave(nodo.segmento):
                der = nodo.segmento
                nodo = nodo.izq
            else:
                izq = nodo.segmento
                nodo = nodo.der

        return izq, der

    # ── en_orden iterativo (stack explícito) ──────────────────────────────────

    def en_orden(self):
        """Devuelve lista de segmentos ordenados por X ascendente en sweep_y."""
        resultado = []
        stack = []
        actual = self._raiz
        while actual or stack:
            while actual:
                stack.append(actual)
                actual = actual.izq
            actual = stack.pop()
            resultado.append(actual.segmento)
            actual = actual.der
        return resultado

    # ── intercambiar ──────────────────────────────────────────────────────────

    def intercambiar(self, seg_a, seg_b):
        self.eliminar(seg_a)
        self.eliminar(seg_b)
        self.insertar(seg_a)
        self.insertar(seg_b)

    # ── utilidades ────────────────────────────────────────────────────────────

    def __len__(self):
        # iterativo con stack
        if not self._raiz:
            return 0
        count = 0
        stack = [self._raiz]
        while stack:
            n = stack.pop()
            count += 1
            if n.izq: stack.append(n.izq)
            if n.der: stack.append(n.der)
        return count

    def __str__(self):
        segs = self.en_orden()
        lineas = [f"  [{i}] {s.p1} -> {s.p2}  (x@sweep={self._x_en_sweep(s):.4f})"
                  for i, s in enumerate(segs)]
        return "ArbolBarrido [\n" + "\n".join(lineas) + "\n]"