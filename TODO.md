# TODO LIST

*(Faites CTRL K+V (je crois) ou CTRL+SHIFT+V pour voir le rendu Markdown sur vsCode)*

Voici une liste non-exhaustive. Je propose, n'hésitez pas à compléter/modfier. L'ordre n'est pas encore décidé. Même si je pense que l'ordre comme tel est plutôt naturel.

Les graphiques ne sont pas la priorité. Ceux-là dépendront des choses à mettre en valeur selon le rapport (format, données, légendes, etc).

## Contexte

Cette section est pour être sûr que l'on parle de la même chose. On a un modèle qui est un **processus de Markov** issu d'un modèle limite *(j'y reviendrai dans quelques semaines pour une l'interprétation d'un résultat, oublions ce "limite")* sur un espace discret fini : $\llbracket 1 ; n \rrbracket$ avec $n$ le nombre de **lignées** ($\approx 50$ pour les simulations numériques).  
Le squelette de ce modèle c'est une chaîne de Markov allant de $n$ à $1$ par pas de $k \in \llbracket 1 ; n-1 \rrbracket$ avec probabilité de $1$.

*La plupart du temps les objets auxquels on s'intéresse dépendent de $n$, à voir si on le rajoute en indice à chaque fois ou si on fait l'impasse et on le laisse implicite.*

Ces lignées peuvent **coalescer** ($k\geq 2$ lignées peuvent fusionner en une seule lignée) avec un taux de coalescence $\lambda_{b,k}$. Cette intensité est définie par une **mesure** $\Lambda$ sur $[0,1]$ Tout tourne autour de celle-là.

Ainsi les problématiques possibles tourneront autour du   
+ **TMRCA** *(Time to the Most Recent Common Ancestor)* (le temps de passer de $n$ lignées à $1$). C'est une variable aléatoire : T_{MRCA}.
+ **Nombre de lignées** à l'instant $t$ noté $N_t$ (on a donc $N_0 = n$ et $N_\infty = 1$).
+ **Nombre de coalescences** jusqu'à l'instant $t$ noté $C_t$ (on a $C_0=0$ et $C_\infty \leq n-1$).
+ **Nombre de k-coalescences** à l'instant $t$ noté $C_t^k$ représente le nombre de coalescences où $k$ lignées ont fusionné en une seule jusqu'à l'instant $t$. On a donc $C_0^k = 0$ et $C_\infty^k \leq \lfloor (n-1)/k \rfloor$ et $\sum_{k=2}^n C_t^k = C_t$
+ **Temps passé** dans chaque état (nombre de lignées) avant la coalescence suivante. On note ça $T_k$ le temps passé dans l'état avec $k$ lignées. (Notation du TD).

Il y a encore à définir mais le gros est là.

## Problématiques



### 2. Mesures proches - stabilité quantitative du TMRCA et de $\mathbb E(T_{MRCA})$

+ Introduction de métriques sur l'espace des mesures : **Wasserstein**, **Total Variation (TV)**, **Hellinger**.
+ Trouver/Ecrire un théorème pour dire
> Soit $\Lambda$ et $\Lambda'$ deux mesures sur $[0,1]$ telles que $d(\Lambda, \Lambda') \leq \epsilon$ pour une métrique $d$ donnée. Alors on a
> $$|\mathbb E_{\Lambda} (T_{MRCA}) - \mathbb E_{\Lambda'} (T_{MRCA})| \leq f(\epsilon)$$
> Ou encore (plus dure je pense) [J'ai mis TV mais à voir ce qu'on mettra comme métrique]
> $$d_{TV} (\mathcal L_{\Lambda} (T_{MRCA}), \mathcal L_{\Lambda'} (T_{MRCA})) \leq f(\epsilon)$$
+ Trouver la fonction $f$ la plus "serrée" possible.
+ Comparer notre estimation numérique avec le théorème et voir à quel point on pourrait améliorer. Tester avec loi $Beta$ et peut-être d'autres.

Cela permet de répondre à cette partie du rapport :  
**``Les différences qualitatives en fonction du choix de la mesure.``**  
**``Comparer les temps au MRCA pour différents $\Lambda$ ``**

### 3. Loi de $N_t$ 

+ Déterminer la loi de $N_t$ pour une mesure $\Lambda$ donnée.

Apparemment ça serait lié à la vitesse de convergence. La densité de $N_t$ nous informerait à propos de la fonction de surivie (i.e queue de distribution, i.e. $\mathbb P(T_{MRCA}>t)$)

### 4. Etudes des événements de coalescence.

+ Déterminer la loi de $C_t$ pour une mesure $\Lambda$ donnée. *(ou alors sa moyenne/variance ? si trop compliqué)*
+ Déterminer (ou au moins empiriquement) la distribution des tailles des coalescences : déterminer la loi de $C_t^k$.
*Par exemple : dans Kingman on a que des coalescences par paires donc $C_t^2 = C_t,\ C_t^k = 0$ pour $k\geq3$.*

La deuxième question est bien plus difficile. *(ChatGPT a une formule apparamment pour $\mathbb E(C_{TMRCA}^k)$* )


### Autres à ajouter / compléter 

Nous devons parler afin de vérifier toutes les consignes de 
+ ``vitesse de coalescence``, je pense que ça revient à étudier **TMRCA** mais peut-être que je me trompe. 

+ ``distribution des temps de coalescence`` : ça peut être le temps entre deux coalescences successives ? Mais dans ce cas c'est juste l'intensité pour toutes les possibilités de coalescences : $Exp(\sum_{k=2}^n \binom{n}{k}\lambda_{n,k})$ donc je ne comprends pas trop ce qui est demandé. Question simple et juste à vérifier empiriquement ?



### Remarques

+ Dans la plupart de la littérature on s'intéresse à des mesures à densités (donc sans atomes). Donc à voir les disjonctions de cas que ça impliquera. Selon la pénibilité on les considerera plus ou moins.

+ A propos du modèle "limite". Le $T_{MRCA}$ est en unité "coalescente". Pour obtenir des unités "réelles" il faut multiplier par une constante dépendant de la population "effective". Or nous sommes dans un modèle limite donc on devra considérer, approximer, que notre nombre de lignées observées $n \ll N_e$ et où on fixera $N_e$ "grand" afin de considérer comme si c'étiat "limite" : $N_e = +\infty$.


# Liens 

+ Papier regardant le TMRCA pour n asymptotique : https://ar5iv.labs.arxiv.org/html/1712.07553
+ Une compilation de résultat probabiliste : https://www.i2m.univ-amu.fr/perso/etienne.pardoux/_media/moehle.pdf
