# MSA_OR-tools
Projet pour le cours de graphes, la complexité et les combinatoires.  

Ce travail a été réalisé dans le cadre du cours graphe complexité et combinatoire du master Intelligence Artificielle de l'université Claude Bernard Lyon 1.  

L'objectif était d'implémenter un modèle en programmation par contrainte d'un problème combinatoire NP-Complet. 

Nous avons choisi de traiter le problème de l'alignement multiple de séquences biologiques (protéiques ou nucléiques) avec le solveur de contraintes OR-Tools de Google. Notre modèle est inspiré du travail de Ignasi Andreu Godall (https://repositori.upf.edu/bitstream/handle/10230/43933/BDBI2019_Andreu_Integer.pdf?sequence=1&isAllowed=y)  

Ici, nous présentons un ensemble de scripts python:
- ```model.py``` : Contient les définitions des modèles et des fonctions utilisées dans les scripts de calculs
- ```MSA.py``` : Réalise les calculs permettant d'avoir les différents alignements en alignements multiples
- ```LCS.py``` : Réalise les calculs permettant d'avoir les plus sous-séquences communes
- ```brute_force.py``` : Calcule les scores de tous les alignements possibles sur 3 séquences et donne l'optimal


---

Project for master's course about Graphs, Complexity, Combinatorials.  

The goal of this work was to implement a model of NP-complete problem into constraint programming or linear programming
