# binary_trees_and_tanglegrams
Sagemath programs related to rooted binary trees and tanglegrams

This package contains Sagemath code to calculate:

1- "minimal_universal_tree_calculator": the list of k-universal rooted binary trees, namely, the rooted binary trees of minimum size, containing every rooted binary tree of size k as an induced subtree.

2- "maximum_number_of_non_isomorphic_subtrees": the number of non-isomorphic induced subtrees of a fixed binary tree. 

3- "tanglegram_isomorphism": Determines whether two tanglegrams are isomorphic, but the algorithm is not shown to work for every size. (works in polynomial time on the number of leaves/vertices.)

4- "tanglegram_multideck_calculator": the k-multideck of a tanglegram T with n leaves, namely, the multiset of tanglegrams "induced" by (n-1)-subsets of the leaves of the left tree, and the corresponding leaves on the right tree, given by the matching. (works in exponential time on the number of leaves/vertices.)

The repository contains the first three programs just mentioned. 

These and related concepts are described in the paper:

Decks of rooted binary trees
Ann Clifton, Eva Czabarka, Audace Dossou-Olory, Kevin Liu, Sarah Loeb, Utku Okur, Laszlo Szekely, Kristina Wicke
https://doi.org/10.48550/arXiv.2311.02255 
(https://arxiv.org/abs/2311.02255)

Abstract: We consider extremal problems related to decks and multidecks of rooted binary trees (a.k.a. rooted phylogenetic tree shapes). Here, the deck (resp. multideck) of a tree T refers to the set (resp. multiset) of leaf induced binary subtrees of T. On the one hand, we consider the reconstruction of trees from their (multi)decks. We give lower and upper bounds on the minimum (multi)deck size required to uniquely encode a rooted binary tree on n leaves. On the other hand, we consider problems related to deck cardinalities. In particular, we characterize trees with minimum-size as well as maximum-size decks. Finally, we present some exhaustive computations for k-universal trees, i.e., rooted binary trees that contain all k-leaf rooted binary trees as induced subtrees.

Relevant Definitions:

1-  A rooted tree is a tree $T$ where a special vertex, denoted by $r_T$, is identified as the root of $T$.
2-  A rooted binary tree is a rooted tree $T$ where every vertex has either no children (i.e., is a leaf ) or two children. 
3-  Let T be a rooted binary tree. For a set $S \subseteq L(T)$, the rooted binary subtree $T[S]$ induced by $S$ is obtained by first taking the minimal connected subgraph $T'$, of $T$ containing $S$, letting $r_{T[S]}$ be the closest vertex of $T'$ to $r_T$, and then suppressing all other degree-2 vertices in $T′$. 

The code is written in Sagemath 10.4, which is available in WSL2 (Windows Subsystem for Linux) in Windows 11.

Download Sagemath: https://www.sagemath.org/

Download WSL2: https://learn.microsoft.com/en-us/windows/wsl/install
