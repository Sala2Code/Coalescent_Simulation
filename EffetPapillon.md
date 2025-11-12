# Effet papillon

Le but est de trouver un théorème stipulant que pour 2 mesures "proches" (de manière quantifiée),
alors l'écart entre leur moyenne de $TMRCA$ est borné finement. Puis, d'appliquer le théorème
sur un cas pratique (loi Beta) et de comparer avec nos simulations numériques l'estimation faite
et le réel écart.

## Rappel

On a montré dans l'introduction de mon rapport,

Partant de $b$ lignées, le taux total auquel une coalescence se produit,
le taux de sortie de l'état $b$, est donné par

$$\begin{equation}
    \lambda_b = \sum_{k=2}^{b} \binom{b}{k} \lambda_{b,k} = \int_0^1 S_b(x) \,\Lambda(dx),
    \quad 
    % \text{avec} \quad 
    S_b(x) \coloneq \sum_{k=2}^{b} \binom{b}{k} x^{k-2}(1-x)^{b-k}
    = \frac{1-(1-x)^b - b x (1-x)^{b-1}}{x^2}
\end{equation}$$
    
Chaque événement de coalescence correspond alors à
une transition de $b$ vers $b-k+1$ lignées avec probabilité 
$$\forall b\geq k \geq 2,
\quad p_{b,k} \coloneq \frac{\binom{b}{k} \lambda_{b,k}}{\lambda_b}
$$

## Concepts

J'y connais rien entre métrique sur les espaces de mesures donc je vais
de me renseigner.
Nous souhaitons comparer des mesures de **probabilité** sur $[0,1]$  
(afin de savoir l'écart résultant des caractéristiques des processus de
coalescence selon 2 mesures proches).
donc une métrique adaptée serait la Total Variation (TV).
*(TODO :Expliquer/Motiver pour celle-là et pas une autre : Hellinger, Wasserstein)*
Elle est lié à KL et Hellinger *(https://fr.wikipedia.org/wiki/Distance_en_variation_totale_(probabilit%C3%A9s))*

Nous travaillons sur $E$, l'espace des mesures de probabilité sur $[0,1]$.
Afin de considérer des mesures proches, nous devons équiper cet espace d'une métrique.
Avec l'analogie de $L_\infty$ sur les fonctionnelles, on définit une métrique associant le plus grand écart
de mesure possible de $E$.

**Definition** : Nous appelons distance en variation totale (TV) entre deux mesures de probabilité $\Lambda_1, \Lambda_2 \in E$,
$$
d_{TV} (\Lambda_1, \Lambda_2) \coloneq \sup_{A\subset [0,1]} |\Lambda_1(A) - \Lambda_2(A)|
$$
**Proposition** : $(E,d_{TV})$ est un espace métrique.  
**Preuve** :  Toutes les propriétes viennent naturellement de $|\cdot|$. Soit $\Lambda_1, \Lambda_2, \Lambda_3 \in E$,
+ $d_{TV}(\Lambda_1, \Lambda_2) = 0 \iff \forall A \subset [0,1], \ \Lambda_1(A) = \Lambda_2(A) \iff \Lambda_1 = \Lambda_2$
+ $d_{TV}(\Lambda_1, \Lambda_2) = d_{TV}(\Lambda_2, \Lambda_1)$
+  $$ d_{TV}(\Lambda_1, \Lambda_3) = \sup_{A\subset [0,1]} |\Lambda_1(A) - \Lambda_2(A) + \Lambda_2(A) - \Lambda_3(A)| \leq d_{TV}(\Lambda_1, \Lambda_2) + d_{TV}(\Lambda_2, \Lambda_3)
$$

Un rappel *(motivier ce rappel, faire un petit blabla)*

**Definition** : Soit $\Lambda, \mu \in E$. On dit que $\Lambda$ est absolument continue par rapport à $\mu$, notée $\lambda \ll \mu$, si pour tout borélien $A \subset [0,1]$ $\mu(A) = 0 \implies \Lambda(A) = 0$.



**Théorème (Radon-Nikodym)** : Soit $\Lambda, \mu \in E$. Si $\Lambda \ll \mu$, alors $\Lambda$ possède une 
densité par rapport à $\mu$.

Cette notion permet d'exprimer la distance de mesures à partir de leurs densités.
**Theoreme** Soit $\Lambda_1, \Lambda_2 \in E$ de densité $u$ et $v$  par rapport $\eta$, une mesure finie telle que $\Lambda_1, \Lambda_2 \ll \eta$. Alors
$$
d_{TV} (\Lambda_1, \Lambda_2) = \frac{1}{2}\|u-v \|_{L_1(\eta)} = \frac{1}{2} \int_{0}^1 |u- v| \, d\eta
$$
**Preuve** Soit $B \subset [0,1]$ un borélien. On a
$$
| \Lambda_1(B)-\Lambda_2(B) | = |\int_B (u-v) \, d\eta | \leq  \int_B |u-v| \, d\eta
$$ 
De plus en utilisant le fait qu'on ait des mesures de probabilité, on obtient de même,
$$
| \Lambda_1(B)-\Lambda_2(B) | = | 1-\Lambda_1(B^C)- (1-\Lambda_2(B^C)) |  = | \Lambda_2(B^C)-\Lambda_1(B^C) |  \leq  \int_{B^C} |u-v| \, d\eta
$$
Donc en sommant et en passant au supremum,
$$
\sup_{B \subset [0,1]} 2 | \Lambda_1(B)-\Lambda_2(B) | \leq \sup_{B \subset [0,1]} \int_B |u-v| \, d\eta + \int_{B^C} |u-v| \, d\eta = \int_0^1 |u-v| \, d\eta
$$

Considérons $A\coloneq \{x \in [0,1],\ u(x)>v(x) \}$ un borélien de $[0,1]$. On a sur $A$ et $A^C$ que $u-v$ est de même signe donc, 
$$
2 | \Lambda_1(A)-\Lambda_2(A) |
= | \int_A (u-v) \, d\eta | + | \int_{A^C} (u-v) \, d\eta |
= \int_A |u-v| \, d\eta + \int_{A^C} |u-v| \, d\eta
= \int_0^1 |u-v| \, d\eta
$$

Le supremum est atteint en $A$. D'où le résultat.