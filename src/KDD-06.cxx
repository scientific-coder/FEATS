/*
Programme pour obtenir les résultats pour l'article à proposer à KDD 06:
"Symbolic Representations of Time Series : a generic framework and evaluations."


Pour les représentations 

Segmentation classique avec modèle linéaire d'ordre 0
Segmentation classique avec modèle linéaire d'ordre 1
SAX
SBSR-L0
CBSR

Donner les infos suivantes
(temps de calcul)
SSE par rapport à une série reconstruite

pour les repr symboliques: SSE par rapport à une "bande" reconstruite pour les symboles extrêmes de SAX, fixer leur aire (aire moyenne des autres symboles de l'alphabet: réalisable a priori, sinon prendre l'aire qui permet de contenir les points extrêmes comme SBSR-L0 ?)

pour les repr symboliques: aire

taille de la repr (taille effective, pour les numéros de symboles, utiliser 8=sizeof(char) plutôt que log2(K), pour les indices sizeof(int) plutôt que log2(N))

pour les repr symb: taille compressée pour tenir compte de la répétition de séquences



*/


/*
charger un ensemble de séries temporelles
transposer celui-ci si necessaire
effectuer un changement de repr et afficher les infos correspondantes pour un tuple de paramètres donné
(nb symbs min, max, mod,nb segs min, max mod)
Pour les segmentations, utilisation d'algo top down pour des raisons de rapidité
*/


/*
faire une nouvelle segmentation equal_length
faire un nouveau model SAX

gestion de l'aire des segments ???
*/
