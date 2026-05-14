import matplotlib.pyplot as plt

def graficar(vertices, aristas):

    for v in vertices.values():
        plt.scatter(v.x, v.y)
        plt.text(v.x, v.y, v.nombre)

    for a in aristas.values():
        origen = vertices[a.origen]
        destino = vertices[aristas[a.antiarista].origen]

        x = [origen.x, destino.x]
        y = [origen.y, destino.y]

        plt.plot(x, y)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()

def viewer(vertices, aristas, caras):

    fig, ax = plt.subplots(figsize=(10, 10))

    # --- aristas ---
    dibujadas = set()
    for a in aristas.values():
        if a.nombre in dibujadas:
            continue
        p1 = a.origen.pt
        p2 = a.antiarista.origen.pt
        ax.plot([p1.x, p2.x], [p1.y, p2.y], color='black', linewidth=1, zorder=3)
        dibujadas.add(a.nombre)
        dibujadas.add(a.antiarista.nombre)

    # --- vertices ---
    for v in vertices.values():
        ax.plot(v.pt.x, v.pt.y, 'ko', markersize=3, zorder=4)


    plt.tight_layout()
    plt.show()