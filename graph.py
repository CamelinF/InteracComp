import matplotlib.pyplot as plt
import igraph as ig
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from main import get_Jaccard

def getLinkFromDico(dictCouv,seuilJacc):
    linksBetweenMotifs=[]
    for i in dictCouv:
        for j in dictCouv:
            if i > j:
                jaccInd=get_Jaccard(dictCouv[i][1],dictCouv[j][1])
                if jaccInd > seuilJacc:
                    if len(dictCouv[j][0]) > len(dictCouv[i][0]):
                        linksBetweenMotifs.append([i,j])
                    else:
                        linksBetweenMotifs.append([j,i])

    return linksBetweenMotifs

def getJaccardLinkFromDico(dictCouv):
    listEdges=[]
    listJaccard=[]
    for i in range(len(CTD)):
        for j in range(len(CTD)):
            if i > j:
                jaccInd=get_Jaccard(dictCouv[i][1],dictCouv[j][1])
                listEdges.append([i,j])
                listJaccard.append(jaccInd)
    return listEdges,listJaccard

def getIntersecLinkFromPattern(CTD):
    listEdges=[]
    listWeigth=[]
    for i in range(len(CTD)):
        for j in range(i+1,len(CTD)):
            m1=set(CTD[i][2])
            m2=set(CTD[j][2])
            inter=m1.intersection(m2)
            listEdges.append([i,j])
            listWeigth.append(len(inter))
    return listEdges,listWeigth


def getLinkFromPattern(CTD, support):
    Links=[]
    for i in range(len(CTD)):
        for j in range(i+1,len(CTD)):
                m1=set(CTD[i][2])
                m2=set(CTD[j][2])
                inter=m1.intersection(m2)
                if inter ==m1:
                    Links.append([j,i])
                elif inter==m2:
                    Links.append([i,j])
                else:
                    if len(inter)>(len(m1)+len(m2))/2: #:support:
                        Links.append([i,j])
    
    return Links


def makeGraph(graph,pr_scores):

       
    # Enlève les noeuds isolés
    # non_isolated_vertices = [v.index for v in graph.vs if graph.degree(v.index) > 0]
    # print(f" on affiche seulements {len(non_isolated_vertices)} sommets")
    # graph = graph.subgraph(non_isolated_vertices)

    # Normalisation manuelle entre 0 et 1 pour l'échelle de couleur
    norm = mcolors.Normalize(vmin=min(pr_scores), vmax=max(pr_scores))
    cmap = cm.get_cmap("coolwarm")  # ou "viridis", "coolwarm", etc.
    colors = [cmap(norm(score)) for score in pr_scores]

    fig, ax = plt.subplots(figsize=(20, 20))
    plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
    layout = graph.layout("fr")

    # Étirement pour espacer les sommets
    scale_factor = 20.0  # Ajuste selon tes besoins
    for coord in layout:
        coord[0] *= scale_factor
        coord[1] *= scale_factor
    vertex_labels = [str(v.index) for v in graph.vs]
    # Dessin du graphe
    ig.plot(
        graph,
        target=ax,
        layout=layout,
        vertex_label=vertex_labels,
        margin=10,
        vertex_color=colors,
        vertex_label_size=8,
        vertex_size=30,
        edge_width=0.7,
        edge_arrow_size=3
    )

    # Ajout de la colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Hack pour éviter un warning
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("PageRank")
    return fig, ax



#  palette = ig.ClusterColoringPalette(len(clusters))
#     colors = [palette[c] for c in g.vs["cluster"]]
#     vertex_labels = [str(v.index) for v in g.vs]
#     # Affichage avec matplotlib
#     layout = g.layout("fr")
#     fig, ax = plt.subplots(figsize=(8, 8))
#     ig.plot(
#         g,
#         target=ax,
#         layout=layout,
#         vertex_label=vertex_labels,
#         vertex_color=colors,
#         vertex_size=30,
#         edge_arrow_size=0.7,
#         margin=30
#     )
#     plt.title(f"Clustering (Louvain) – {len(clusters)} communautés")
#     plt.show()