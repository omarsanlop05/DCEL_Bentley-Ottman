# -*- coding: utf-8 -*-

import pygame


# ============================================================
# PALETA
# ============================================================

_PALETA = [
    (173, 216, 230),
    (144, 238, 144),
    (255, 255, 153),
    (255, 179, 128),
    (216, 191, 216),
    (255, 182, 193),
    (152, 251, 152),
    (135, 206, 250),
    (240, 230, 140),
    (221, 160, 221),
    (175, 238, 238),
    (255, 218, 185),
]


# ============================================================
# GEOMETRÍA
# ============================================================

def _en_poligono(px, py, poli):

    dentro = False

    j = len(poli) - 1

    for i in range(len(poli)):

        xi, yi = poli[i]
        xj, yj = poli[j]

        intersecta = (
            ((yi > py) != (yj > py))
            and
            (
                px <
                (xj - xi) * (py - yi)
                / (yj - yi + 1e-12)
                + xi
            )
        )

        if intersecta:
            dentro = not dentro

        j = i

    return dentro


def _area_poligono(puntos):

    area = 0

    for i in range(len(puntos)):

        x1, y1 = puntos[i]
        x2, y2 = puntos[(i + 1) % len(puntos)]

        area += (x1 * y2 - x2 * y1)

    return abs(area) / 2.0


def _centroide(puntos):

    cx = sum(p[0] for p in puntos) / len(puntos)
    cy = sum(p[1] for p in puntos) / len(puntos)

    return cx, cy


# ============================================================
# VISOR
# ============================================================

class VisorDCEL:

    def __init__(
            self,
            vertices_g,
            aristas_g,
            caras_g,
            ancho=1100,
            alto=900,
            titulo="Visor DCEL"
    ):

        self.vertices_g = vertices_g
        self.aristas_g = aristas_g
        self.caras_g = caras_g

        self.ancho = ancho
        self.alto = alto
        self.titulo = titulo

        # ----------------------------------------------------
        # cámara
        # ----------------------------------------------------

        self._escala = 1.0
        self._offset_x = 0
        self._offset_y = 0

        # ----------------------------------------------------
        # interacción
        # ----------------------------------------------------

        self._seleccionadas = {}

        self._color_idx = 0

        self._pan_activo = False

        self._pan_inicio = (0, 0)

        self._pan_offset_0 = (0, 0)

        # ----------------------------------------------------
        # geometría
        # ----------------------------------------------------

        self._caras_data = []

        self._aristas_segs = []

        self._extraer_geometria()

        self._calcular_holes_geometricos()

        self._ajustar_vista()

    # ========================================================
    # EXTRACCIÓN
    # ========================================================

    def _extraer_geometria(self):

        # ----------------------------------------------------
        # CARAS
        # ----------------------------------------------------

        for nombre, cara in self.caras_g.items():

            if nombre == "f_infinita":
                continue

            if cara.aristas_externas is None:
                continue

            puntos = []

            actual = cara.aristas_externas

            visitadas = set()

            while actual and actual.nombre not in visitadas:

                visitadas.add(actual.nombre)

                puntos.append((
                    actual.origen.pt.x,
                    actual.origen.pt.y
                ))

                actual = actual.siguiente

            if len(puntos) < 3:
                continue

            self._caras_data.append({

                "nombre": nombre,

                "puntos": puntos,

                "centroide": _centroide(puntos),

                "area": _area_poligono(puntos),

                "holes": [],
            })

        # ----------------------------------------------------
        # ordenar pequeñas → grandes
        # ----------------------------------------------------

        self._caras_data.sort(
            key=lambda c: c["area"]
        )

        # ----------------------------------------------------
        # ARISTAS
        # ----------------------------------------------------

        vistas = set()

        for nombre, a in self.aristas_g.items():

            if nombre in vistas:
                continue

            p1 = (
                a.origen.pt.x,
                a.origen.pt.y
            )

            p2 = (
                a.antiarista.origen.pt.x,
                a.antiarista.origen.pt.y
            )

            self._aristas_segs.append((p1, p2))

            vistas.add(nombre)
            vistas.add(a.antiarista.nombre)

    # ========================================================
    # HOLES GEOMÉTRICOS AUTOMÁTICOS
    # ========================================================

    def _calcular_holes_geometricos(self):

        n = len(self._caras_data)

        for i in range(n):

            cara_grande = self._caras_data[i]

            poli_grande = cara_grande["puntos"]

            for j in range(i):

                cara_pequena = self._caras_data[j]

                cx, cy = cara_pequena["centroide"]

                if _en_poligono(cx, cy, poli_grande):

                    cara_grande["holes"].append(
                        cara_pequena["puntos"]
                    )

    # ========================================================
    # VISTA
    # ========================================================

    def _ajustar_vista(self):

        xs = [
            p[0]
            for c in self._caras_data
            for p in c["puntos"]
        ]

        ys = [
            p[1]
            for c in self._caras_data
            for p in c["puntos"]
        ]

        min_x = min(xs)
        max_x = max(xs)

        min_y = min(ys)
        max_y = max(ys)

        margen = 50

        rango_x = max_x - min_x or 1
        rango_y = max_y - min_y or 1

        self._escala = min(
            (self.ancho - margen * 2) / rango_x,
            (self.alto - margen * 2) / rango_y
        )

        cx = (min_x + max_x) / 2
        cy = (min_y + max_y) / 2

        self._offset_x = self.ancho / 2 - cx * self._escala
        self._offset_y = self.alto / 2 + cy * self._escala

    # ========================================================
    # TRANSFORMACIONES
    # ========================================================

    def _m2p(self, x, y):

        return (
            int(x * self._escala + self._offset_x),
            int(-y * self._escala + self._offset_y)
        )

    def _poli_px(self, puntos):

        return [
            self._m2p(x, y)
            for x, y in puntos
        ]

    # ========================================================
    # HIT TEST
    # ========================================================

    def _cara_en(self, px, py):

        candidatas = []

        for cara in self._caras_data:

            poli = self._poli_px(cara["puntos"])

            if _en_poligono(px, py, poli):

                candidatas.append(cara)

        if not candidatas:
            return None

        return min(
            candidatas,
            key=lambda c: c["area"]
        )

    # ========================================================
    # DIBUJO
    # ========================================================

    def _dibujar(self, surface, fuente_nombre, fuente_ayuda):

        surface.fill((255, 255, 255))

        # ----------------------------------------------------
        # CARAS
        # ----------------------------------------------------

        for cara in self._caras_data:

            if cara["nombre"] not in self._seleccionadas:
                continue

            color = self._seleccionadas[cara["nombre"]]

            temp = pygame.Surface(
                (self.ancho, self.alto),
                pygame.SRCALPHA
            )

            # exterior
            pygame.draw.polygon(
                temp,
                (*color, 255),
                self._poli_px(cara["puntos"])
            )

            # holes geométricos
            for hole in cara["holes"]:

                pygame.draw.polygon(
                    temp,
                    (0, 0, 0, 0),
                    self._poli_px(hole)
                )

            surface.blit(temp, (0, 0))

        # ----------------------------------------------------
        # ARISTAS
        # ----------------------------------------------------

        for (x1, y1), (x2, y2) in self._aristas_segs:

            pygame.draw.line(
                surface,
                (30, 30, 30),
                self._m2p(x1, y1),
                self._m2p(x2, y2),
                1
            )

        # ----------------------------------------------------
        # NOMBRES
        # ----------------------------------------------------

        for cara in self._caras_data:

            if cara["nombre"] not in self._seleccionadas:
                continue

            cx, cy = self._m2p(*cara["centroide"])

            txt = fuente_nombre.render(
                cara["nombre"],
                True,
                (20, 20, 20)
            )

            rect = txt.get_rect(center=(cx, cy))

            surface.blit(txt, rect)

        # ----------------------------------------------------
        # HUD
        # ----------------------------------------------------

        ayuda = [
            "Click izq : seleccionar",
            "Rueda     : zoom",
            "Click der : pan",
            "R : reset",
        ]

        y0 = self.alto - 90

        for linea in ayuda:

            surf = fuente_ayuda.render(
                linea,
                True,
                (100, 100, 100)
            )

            surface.blit(surf, (10, y0))

            y0 += 18

    # ========================================================
    # LOOP
    # ========================================================

    def ejecutar(self):

        pygame.init()

        screen = pygame.display.set_mode(
            (self.ancho, self.alto)
        )

        pygame.display.set_caption(self.titulo)

        clock = pygame.time.Clock()

        fuente_nombre = pygame.font.SysFont(
            "Arial",
            13,
            bold=True
        )

        fuente_ayuda = pygame.font.SysFont(
            "Arial",
            12
        )

        corriendo = True

        while corriendo:

            for ev in pygame.event.get():

                # ------------------------------------------------
                # salir
                # ------------------------------------------------

                if ev.type == pygame.QUIT:
                    corriendo = False

                # ------------------------------------------------
                # teclado
                # ------------------------------------------------

                elif ev.type == pygame.KEYDOWN:

                    if ev.key == pygame.K_ESCAPE:
                        corriendo = False

                    elif ev.key == pygame.K_r:

                        self._seleccionadas.clear()

                        self._color_idx = 0

                # ------------------------------------------------
                # zoom
                # ------------------------------------------------

                elif ev.type == pygame.MOUSEWHEEL:

                    mx, my = pygame.mouse.get_pos()

                    factor = 1.15 if ev.y > 0 else 1 / 1.15

                    self._offset_x = (
                        mx - (mx - self._offset_x) * factor
                    )

                    self._offset_y = (
                        my - (my - self._offset_y) * factor
                    )

                    self._escala *= factor

                # ------------------------------------------------
                # mouse down
                # ------------------------------------------------

                elif ev.type == pygame.MOUSEBUTTONDOWN:

                    # PAN

                    if ev.button == 3:

                        self._pan_activo = True

                        self._pan_inicio = ev.pos

                        self._pan_offset_0 = (
                            self._offset_x,
                            self._offset_y
                        )

                    # SELECCIÓN

                    elif ev.button == 1:

                        cara = self._cara_en(*ev.pos)

                        if cara:

                            nombre = cara["nombre"]

                            if nombre in self._seleccionadas:

                                del self._seleccionadas[nombre]

                            else:

                                self._seleccionadas[nombre] = (
                                    _PALETA[
                                        self._color_idx
                                        % len(_PALETA)
                                    ]
                                )

                                self._color_idx += 1

                # ------------------------------------------------
                # mouse up
                # ------------------------------------------------

                elif ev.type == pygame.MOUSEBUTTONUP:

                    if ev.button == 3:
                        self._pan_activo = False

                # ------------------------------------------------
                # pan move
                # ------------------------------------------------

                elif ev.type == pygame.MOUSEMOTION:

                    if self._pan_activo:

                        dx = ev.pos[0] - self._pan_inicio[0]
                        dy = ev.pos[1] - self._pan_inicio[1]

                        self._offset_x = (
                            self._pan_offset_0[0] + dx
                        )

                        self._offset_y = (
                            self._pan_offset_0[1] + dy
                        )

            self._dibujar(
                screen,
                fuente_nombre,
                fuente_ayuda
            )

            pygame.display.flip()

            clock.tick(60)

        pygame.quit()