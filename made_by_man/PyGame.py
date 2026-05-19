# -*- coding: utf-8 -*-
"""


Controles:
    Clic izquierdo           → selecciona / deselecciona una cara
    Rueda del mouse          → zoom
    Clic derecho + arrastrar → mover el mapa (pan)
    R                        → reinicia todos los colores
    ESC / cerrar             → salir
"""

import pygame
import sys


# ─────────────────────────────────────────────────────────────────────────────
# Paleta pastel
# ─────────────────────────────────────────────────────────────────────────────
_PALETA = [
    (173, 216, 230), (144, 238, 144), (255, 255, 153), (255, 179, 128),
    (216, 191, 216), (255, 182, 193), (152, 251, 152), (135, 206, 250),
    (240, 230, 140), (221, 160, 221), (175, 238, 238), (255, 218, 185),
    (204, 255, 204), (255, 204, 229), (204, 229, 255), (255, 229, 204),
]


def _en_poligono(px, py, poli):
    """Ray-casting: True si el punto (px,py) está dentro del polígono."""
    dentro = False
    n = len(poli)
    j = n - 1
    for i in range(n):
        xi, yi = poli[i]
        xj, yj = poli[j]
        if (yi > py) != (yj > py):
            ix = (xj - xi) * (py - yi) / (yj - yi + 1e-12) + xi
            if px < ix:
                dentro = not dentro
        j = i
    return dentro


class VisorDCEL:
    """
    Ventana Pygame para explorar la DCEL de forma interactiva.

    Parámetros
    ----------
    vertices_g, aristas_g, caras_g : dicts producidos por main.py
    ancho, alto                    : tamaño de la ventana en píxeles
    titulo                         : texto de la barra de título
    """

    def __init__(self, vertices_g, aristas_g, caras_g,
                 ancho=950, alto=800, titulo="Visor DCEL Interactivo"):

        self.vertices_g = vertices_g
        self.aristas_g  = aristas_g
        self.caras_g    = caras_g
        self.ancho      = ancho
        self.alto       = alto
        self.titulo     = titulo

        # Transformación mundo → pantalla
        self._escala   = 1.0
        self._offset_x = 0.0
        self._offset_y = 0.0

        # Estado interactivo
        self._seleccionadas = {}   # nombre_cara -> (r, g, b)
        self._cara_exterior_nombre = "__EXTERIOR__"
        self._color_idx     = 0
        self._pan_activo    = False
        self._pan_inicio    = (0, 0)
        self._pan_offset_0  = (0.0, 0.0)

        # Geometría pre-calculada en coords mundo
        self._caras_data   = []   # [{ nombre, puntos, centroide }]
        self._aristas_segs = []   # [ ((x1,y1),(x2,y2)), ... ]

        self._extraer_geometria()
        self._ajustar_vista()

    # ─── extracción ──────────────────────────────────────────────────────────

    def _extraer_geometria(self):
        # Caras
        for nombre, cara in self.caras_g.items():
            if nombre == "f_infinita" or cara.aristas_externas is None:
                continue
            puntos    = []
            actual    = cara.aristas_externas
            visitadas = set()
            while actual and actual.nombre not in visitadas:
                visitadas.add(actual.nombre)
                puntos.append((actual.origen.pt.x, actual.origen.pt.y))
                actual = actual.siguiente
            if len(puntos) >= 3:
                cx = sum(p[0] for p in puntos) / len(puntos)
                cy = sum(p[1] for p in puntos) / len(puntos)
                self._caras_data.append({
                    "nombre":    nombre,
                    "puntos":    puntos,
                    "centroide": (cx, cy),
                })

        # Aristas (sin duplicar la media-arista gemela)
        vistas = set()
        for nombre, a in self.aristas_g.items():
            if nombre in vistas:
                continue
            p1 = (a.origen.pt.x,           a.origen.pt.y)
            p2 = (a.antiarista.origen.pt.x, a.antiarista.origen.pt.y)
            self._aristas_segs.append((p1, p2))
            vistas.add(nombre)
            vistas.add(a.antiarista.nombre)

    def _ajustar_vista(self):
        """Escala y centra el mapa completo en la ventana."""
        if not self._caras_data:
            return
        todos_x = [p[0] for c in self._caras_data for p in c["puntos"]]
        todos_y = [p[1] for c in self._caras_data for p in c["puntos"]]
        min_x, max_x = min(todos_x), max(todos_x)
        min_y, max_y = min(todos_y), max(todos_y)

        margen   = 50
        rango_x  = max_x - min_x or 1
        rango_y  = max_y - min_y or 1
        self._escala = min(
            (self.ancho - 2 * margen) / rango_x,
            (self.alto  - 2 * margen) / rango_y,
        )
        cx_m = (min_x + max_x) / 2
        cy_m = (min_y + max_y) / 2
        self._offset_x = self.ancho / 2 - cx_m *  self._escala
        self._offset_y = self.alto  / 2 + cy_m *  self._escala   # Y invertido

    # ─── coordenadas ─────────────────────────────────────────────────────────

    def _m2p(self, x, y):
        """Mundo → píxeles (Y invertido)."""
        return (
            int(x *  self._escala + self._offset_x),
            int(y * -self._escala + self._offset_y),
        )

    def _poli_px(self, puntos_mundo):
        return [self._m2p(x, y) for x, y in puntos_mundo]

    # ─── hit-test ────────────────────────────────────────────────────────────

    def _cara_en(self, px, py):
        """
        Devuelve la cara bajo el cursor.

        Si el punto no cae dentro de ninguna cara interna,
        entonces pertenece a la cara exterior (f0).
        """
        for cara in reversed(self._caras_data):
            if _en_poligono(px, py, self._poli_px(cara["puntos"])):
                return cara

        return {
            "nombre": self._cara_exterior_nombre,
            "puntos": [],
            "centroide": (0, 0),
        }

    # ─── dibujo ──────────────────────────────────────────────────────────────

    def _dibujar(self, surface, fuente_nombre, fuente_ayuda):

        surface.fill((255, 255, 255))

        if self._cara_exterior_nombre in self._seleccionadas:
            color_ext = self._seleccionadas[self._cara_exterior_nombre]
            # 1. Pintamos TODO
            surface.fill(color_ext)
            # 2. "Recortamos" las caras internas
            for cara in self._caras_data:
                poli = self._poli_px(cara["puntos"])
                if len(poli) >= 3:
                    pygame.draw.polygon(surface, (255, 255, 255), poli)


        # 1. Relleno de caras seleccionadas
        for cara in self._caras_data:
            if cara["nombre"] not in self._seleccionadas:
                continue
            poli = self._poli_px(cara["puntos"])
            if len(poli) >= 3:
                pygame.draw.polygon(surface, self._seleccionadas[cara["nombre"]], poli)

        # 2. Todas las aristas en negro
        for (x1, y1), (x2, y2) in self._aristas_segs:
            pygame.draw.line(surface, (30, 30, 30),
                             self._m2p(x1, y1), self._m2p(x2, y2), 1)

        # 3. Nombre de cada cara seleccionada sobre su centroide
        for cara in self._caras_data:
            if cara["nombre"] not in self._seleccionadas:
                continue
            cx, cy = self._m2p(*cara["centroide"])
            txt  = fuente_nombre.render(cara["nombre"], True, (20, 20, 20))
            rect = txt.get_rect(center=(cx, cy))
            # fondo blanco semitransparente
            pad = 3
            bg  = pygame.Surface((rect.width + pad*2, rect.height + pad*2),
                                  pygame.SRCALPHA)
            bg.fill((255, 255, 255, 180))
            surface.blit(bg,  (rect.x - pad, rect.y - pad))
            surface.blit(txt, rect)

        # 4. HUD de ayuda (abajo a la izquierda)
        ayuda = [
            "Clic izq : seleccionar / deseleccionar cara",
            "Rueda    : zoom",
            "Clic der + arrastrar : mover mapa",
            "R : limpiar selección   |   ESC : salir",
        ]
        y0 = self.alto - len(ayuda) * 17 - 6
        for linea in ayuda:
            surf = fuente_ayuda.render(linea, True, (120, 120, 120))
            surface.blit(surf, (8, y0))
            y0 += 17

        # 5. Contador (arriba a la derecha)
        conteo = fuente_ayuda.render(
            f"Seleccionadas: {len(self._seleccionadas)} / {len(self._caras_data) + 1}",
            True, (60, 60, 60))
        surface.blit(conteo, (self.ancho - conteo.get_width() - 10, 8))

    # ─── bucle ───────────────────────────────────────────────────────────────

    def ejecutar(self):
        pygame.init()
        screen = pygame.display.set_mode((self.ancho, self.alto))
        pygame.display.set_caption(self.titulo)
        clock  = pygame.time.Clock()

        fuente_nombre = pygame.font.SysFont("Arial", 13, bold=True)
        fuente_ayuda  = pygame.font.SysFont("Arial", 12)

        corriendo = True
        while corriendo:
            for ev in pygame.event.get():

                if ev.type == pygame.QUIT:
                    corriendo = False

                elif ev.type == pygame.KEYDOWN:
                    if ev.key == pygame.K_ESCAPE:
                        corriendo = False
                    elif ev.key == pygame.K_r:
                        self._seleccionadas.clear()
                        self._color_idx = 0

                elif ev.type == pygame.MOUSEWHEEL:
                    mx, my = pygame.mouse.get_pos()
                    factor = 1.15 if ev.y > 0 else 1 / 1.15
                    # zoom anclado al cursor
                    self._offset_x = mx - (mx - self._offset_x) * factor
                    self._offset_y = my - (my - self._offset_y) * factor
                    self._escala  *= factor

                elif ev.type == pygame.MOUSEBUTTONDOWN:
                    if ev.button == 3:   # clic derecho → inicio del pan
                        self._pan_activo   = True
                        self._pan_inicio   = ev.pos
                        self._pan_offset_0 = (self._offset_x, self._offset_y)

                    elif ev.button == 1:   # clic izquierdo → selección
                        cara = self._cara_en(*ev.pos)
                        if cara:
                            nombre = cara["nombre"]
                            if nombre in self._seleccionadas:
                                del self._seleccionadas[nombre]
                            else:
                                color = _PALETA[self._color_idx % len(_PALETA)]
                                self._seleccionadas[nombre] = color
                                self._color_idx += 1

                elif ev.type == pygame.MOUSEBUTTONUP:
                    if ev.button == 3:
                        self._pan_activo = False

                elif ev.type == pygame.MOUSEMOTION:
                    if self._pan_activo:
                        dx = ev.pos[0] - self._pan_inicio[0]
                        dy = ev.pos[1] - self._pan_inicio[1]
                        self._offset_x = self._pan_offset_0[0] + dx
                        self._offset_y = self._pan_offset_0[1] + dy

            self._dibujar(screen, fuente_nombre, fuente_ayuda)
            pygame.display.flip()
            clock.tick(60)

        pygame.quit()