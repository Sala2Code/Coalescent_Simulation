import numpy as np
import matplotlib.pyplot as plt
import numpy.random as npr

# ==============================================================



def simulate_kingman_tree(n):
    k = n
    t = 0.0

    active = [f"leaf_{i}" for i in range(n)]
    tree = []
    node_id = 0

    while k > 1:
        rate = k * (k - 1) / 2
        dt = (-1/rate)*np.log(1-npr.rand())
        t += dt

        # Choisir 2 lignées au hasard à fusionner
        children = npr.choice(active, size=2, replace=False)
        parent = f"node_{node_id}"
        node_id += 1

        tree.append((parent, list(children), t))

        active = [x for x in active if x not in children] + [parent]
        k -= 1

    root = active[0]
    return tree, root

#####################################################3

def plot_tree_optimized(tree, root, n, title="faut bien mettre un titre"):
    positions = {f"leaf_{i}": (i, 0) for i in range(n)}

    def compute_positions(node):
        for parent, children, t in tree:
            if parent == node:
                xs = []
                for c in children:
                    if c not in positions:
                        compute_positions(c)
                    xs.append(positions[c][0])
                x_parent = np.mean(xs)
                positions[parent] = (x_parent, t)

    compute_positions(root)

    plt.figure(figsize=(11, 6))
    
    for parent, children, t in tree:
        xs = []
        for c in children:
            x_c, y_c = positions[c]
            xs.append(x_c)
            plt.plot([x_c, x_c], [y_c, t], color="b", lw=0.5)
        plt.plot([min(xs), max(xs)], [t, t], color="r", lw=2)

    plt.scatter(np.arange(0,n), np.zeros(n), color="black", s=10)

    plt.title(title, fontsize=13)
    plt.xlabel("Individus initiaux")
    plt.ylabel("Temps de coalescence (vers le passé)")
    plt.gca().invert_yaxis()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


# ==============================================================

def simulate_kingman_tree_nonvisual(n):
    k = n
    t = 0.0

    active = [f"leaf_{i}" for i in range(n)]
    tree = []
    node_id = 0

    while k > 1:
        rate = k * (k - 1) / 2
        dt = (-1/rate)*np.log(1-npr.rand())
        t += dt

        # Choisir 2 lignées au hasard à fusionner
        children = active[0], active[1]
        parent = f"node_{node_id}"
        node_id += 1

        tree.append((parent, list(children), t))

        active = active[2:] + [parent]
        k -= 1

    root = active[0]
    return tree, root


def distribution(n, n_simulations=1000):
    total_times = []

    for _ in range(n_simulations):
        tree, root = simulate_kingman_tree_nonvisual(n)
        total_time = tree[-1][2]  # Temps du dernier événement de coalescence
        total_times.append(total_time)
    
    plt.hist(total_times, bins=np.linspace(np.min(total_times),6,30), density=True)

    return np.array(total_times)

n = 500
#tree, root = simulate_kingman_tree(n)
#plot_tree_optimized(tree, root, n, title=f"Arbre de coalescence de Kingman (n={n})")
times = distribution(n, n_simulations=5000)
print(times)
