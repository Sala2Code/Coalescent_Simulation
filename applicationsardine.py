import numpy as np
import matplotlib.pyplot as plt
import numpy.random as npr
from scipy.special import comb, gammaln

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

#####################################################
def simulate_coalescent_with_pairwise(n, alpha, expansion_rate=0):
    active_lineages = [[i] for i in range(n)]
    pairwise_times = []
    t = 0
    
    while len(active_lineages) > 1:
        b = len(active_lineages)
        

        if alpha >= 1.99: # Kingman
            base_rate = b * (b - 1) / 2.0
            probs = None
        else: # Beta (Sardine)
            rates = compute_rates_beta(b, alpha)
            base_rate = rates.sum()
            probs = rates / base_rate


        # --- Ajustement pour Expansion Démographique ---
        current_pop_factor = np.exp(-expansion_rate * t) 
        
        effective_rate = base_rate * (1.0 / current_pop_factor)

        dt = npr.exponential(1.0 / effective_rate)
        t += dt
        
        # --- Choix des lignées à fusionner ---
        if alpha >= 1.98:
            k = 2
            indices = np.random.choice(b, size=2, replace=False)
        else:
            k = np.random.choice(np.arange(2, b + 1), p=probs)
            indices = np.random.choice(b, size=k, replace=False)
        

        merging_lineages = [active_lineages[i] for i in indices]
        
        # Pour toutes les paires qui fusionnent ici, la distance temporelle est t
        # On optimise en calculant juste le nombre de paires
        sizes = [len(l) for l in merging_lineages]
        total_pairs_in_merge = 0
        
        # Somme des produits des tailles (paires entre groupes distincts)
        # Exemple: Groupe A (taille 2) fusionne avec B (taille 3). 2*3 = 6 paires coalescent.
        for i in range(len(sizes)):
            for j in range(i+1, len(sizes)):
                pairwise_times.extend([t] * (sizes[i] * sizes[j]))
        
        new_lineage = []
        for lin in merging_lineages:
            new_lineage.extend(lin)
            
        next_active = [active_lineages[i] for i in range(b) if i not in indices]
        next_active.append(new_lineage)
        active_lineages = next_active

    return np.array(pairwise_times)

####################################################3
def run_sardine_study():
    N = 100
    N_SIM = 100
    TARGET_DIFF = 10    # On vise 10 différences en moyenne pour étaler le graphe
    
    print(f"Simulation : {N_SIM} arbres de {N} feuilles.")

    times_kingman = []
    for _ in range(N_SIM):
        #expansion=5.0 (Forte croissance necessaire)
        times = simulate_coalescent_with_pairwise(N, alpha=2.0, expansion_rate=5.0)#alpha = 2 ne pose pas de probleme, c est simulé differemment
        times_kingman.extend(times)
    times_kingman = np.array(times_kingman)

    print("simu 1 finie")
    times_sardine = []
    for _ in range(N_SIM):
        times = simulate_coalescent_with_pairwise(N, alpha=1.3, expansion_rate=0)
        times_sardine.extend(times)
    times_sardine = np.array(times_sardine)


    mu_kingman = TARGET_DIFF / ( np.mean(times_kingman))
    mu_sardine = TARGET_DIFF / ( np.mean(times_sardine))
    # Ajout du bruit de Poisson (Mutation process)
    mm_kingman = np.random.poisson(mu_kingman *  times_kingman)
    mm_sardine = np.random.poisson(mu_sardine *  times_sardine)

    plt.figure(figsize=(12, 7))
    bins = np.arange(0, 40, 1)
    plt.hist(mm_kingman, bins=bins, density=True, alpha=0.5, color='blue', label='Kingman + Expansion ', edgecolor='black')
    plt.hist(mm_sardine, bins=bins, density=True, alpha=0.5, color='red', label='Sardine (Beta-Coalescent)', edgecolor='darkred')

    plt.axvline(np.mean(mm_kingman), color='blue', linestyle='--', alpha=0.5, label = "Nombre moyen de mutations")
    
    plt.title("La Preuve par la Discordance (Mismatch Distribution)\nLa 'Vague' (Bleu) vs La Réalité (Rouge)")
    plt.xlabel("Nombre de différences génétiques (Mutations)")
    plt.ylabel("Taux")
    #plt.yscale('log')
    
    plt.legend()
    plt.tight_layout()
    plt.show()

run_sardine_study()