/*
Programme pour obtenir les r�sultats pour l'article � proposer � KDD 06:
"Symbolic Representations of Time Series : a generic framework and evaluations."


Pour les repr�sentations 

Segmentation classique avec mod�le lin�aire d'ordre 0
Segmentation classique avec mod�le lin�aire d'ordre 1
SAX
SBSR-L0
CBSR

Donner les infos suivantes
(temps de calcul)
SSE par rapport � une s�rie reconstruite

pour les repr symboliques: SSE par rapport � une "bande" reconstruite pour les symboles extr�mes de SAX, fixer leur aire (aire moyenne des autres symboles de l'alphabet: r�alisable a priori, sinon prendre l'aire qui permet de contenir les points extr�mes comme SBSR-L0 ?)

pour les repr symboliques: aire

taille de la repr (taille effective, pour les num�ros de symboles, utiliser 8=sizeof(char) plut�t que log2(K), pour les indices sizeof(int) plut�t que log2(N))

pour les repr symb: taille compress�e pour tenir compte de la r�p�tition de s�quences



*/


/*
charger un ensemble de s�ries temporelles
transposer celui-ci si necessaire
effectuer un changement de repr et afficher les infos correspondantes pour un tuple de param�tres donn�
(nb symbs min, max, mod,nb segs min, max mod)
Pour les segmentations, utilisation d'algo top down pour des raisons de rapidit�
*/


/*
faire une nouvelle segmentation equal_length
faire un nouveau model SAX

gestion de l'aire des segments ???
*/
