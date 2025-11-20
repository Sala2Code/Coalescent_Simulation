import numpy as np
import matplotlib.pyplot as plt
import numpy.random as npr
from scipy.special import gamma, comb, gammaln

def compute_rates_beta(b, alpha): #on passe alors en espace log pour eviter les problemes de flottants precition et ça nous permet de prendre des n plus grands
    
    """
    Computes rates for Beta-Coalescent (0 < alpha < 2).
    Formula: Lambda_{b,k} =  (gamma(k-alpha)*gamma(b-k+alpha))/(gamma(b)*gamma(2-alpha)*gamma(alpha))
    """
    k = np.arange(2, b + 1)
    

    lognum = gammaln(k-alpha) + gammaln(b - k + alpha)
    logdenom = gammaln(2 - alpha) + gammaln(alpha) + gammaln(b)
    
    log_rates = lognum - logdenom
    return comb(b,k) * np.exp(log_rates)

def lambda_bk(b,k, alpha): #formule basique premier essai, correcte mais echoue n > 50 car probleme folattants
    return (gamma(k-alpha)*gamma(b-k+alpha))/(gamma(b)*gamma(2-alpha)*gamma(alpha))


def lambda_bk_montecarlo(b, k, alpha, N=10000): #fonctionne pas pour l instant
    x = npr.rand(N)
    values = x**(k-1+alpha) * (1-x)**(b-k+alpha-1)
    return np.mean(values)

def simulate_coalescent(n, alpha):
    active = [f"leaf_{i}" for i in range(n)]
    tree = []
    node_id = 0
    t = 0 
    while len(active) > 1:
        b = len(active)
        
        
        rates = compute_rates_beta(b, alpha)
        total_rate = rates.sum()
        
        # temps jusqu'à la prochaine fusion (exponentiel)
        dt = np.random.exponential(1)/total_rate
        t += dt
        
        probs = rates / total_rate

        k = np.random.choice(range(2,b+1), p=probs)


        children = npr.choice(active, size=k, replace=False)
        parent = f"node_{node_id}"
        node_id += 1

        tree.append((parent, list(children), t))

        chosen_set = set(children)
        active = [x for x in active if x not in chosen_set]
        active.append(parent)

        

    root = active[0]
    return tree, root
print(simulate_coalescent(10, 1.5))


######################################################
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

########################################################
def distribution(n, n_simulations=1000):
    total_times = []

    for _ in range(n_simulations):
        tree, root = simulate_coalescent(n, alpha = 1)
        total_time = tree[-1][2]  # Temps du dernier événement de coalescence
        total_times.append(total_time)
    
    plt.hist(total_times, bins=np.linspace(np.min(total_times),6,30), density=True)

    return np.array(total_times)
n = 1000
tree, root = simulate_coalescent(n, alpha = 1.5)
plot_tree_optimized(tree, root, n, title=f"Arbre de coalescence  beta (n={n})")
plt.show()
#times = distribution(n/2, n_simulations=5000)