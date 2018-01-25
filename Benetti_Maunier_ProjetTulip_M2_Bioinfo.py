# Written by Benetti Julien & Maunier Tristan
# Powered by Python 3.6

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import tlp
import math
from scipy.stats.stats import pearsonr 

# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# The runGraphScript(scriptFile, graph) function can be called to launch
# another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# The main(graph) function must be defined 
# to run the script on the current graph

def main(graph): 
  Locus = graph.getStringProperty("Locus")
  Negative = graph.getBooleanProperty("Negative")
  Positive = graph.getBooleanProperty("Positive")
  tp1_s = graph.getDoubleProperty("tp1 s")
  tp10_s = graph.getDoubleProperty("tp10 s")
  tp11_s = graph.getDoubleProperty("tp11 s")
  tp12_s = graph.getDoubleProperty("tp12 s")
  tp13_s = graph.getDoubleProperty("tp13 s")
  tp14_s = graph.getDoubleProperty("tp14 s")
  tp15_s = graph.getDoubleProperty("tp15 s")
  tp16_s = graph.getDoubleProperty("tp16 s")
  tp17_s = graph.getDoubleProperty("tp17 s")
  tp2_s = graph.getDoubleProperty("tp2 s")
  tp3_s = graph.getDoubleProperty("tp3 s")
  tp4_s = graph.getDoubleProperty("tp4 s")
  tp5_s = graph.getDoubleProperty("tp5 s")
  tp6_s = graph.getDoubleProperty("tp6 s")
  tp7_s = graph.getDoubleProperty("tp7 s")
  tp8_s = graph.getDoubleProperty("tp8 s")
  tp9_s = graph.getDoubleProperty("tp9 s")
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
  viewColor = graph.getColorProperty("viewColor")
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  viewLabel = graph.getStringProperty("viewLabel")
  viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
  viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
  viewLabelColor = graph.getColorProperty("viewLabelColor")
  viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewMetric = graph.getDoubleProperty("viewMetric")
  viewRotation = graph.getDoubleProperty("viewRotation")
  viewSelection = graph.getBooleanProperty("viewSelection")
  viewShape = graph.getIntegerProperty("viewShape")
  viewSize = graph.getSizeProperty("viewSize")
  viewSrcAnchorShape = graph.getIntegerProperty("viewSrcAnchorShape")
  viewSrcAnchorSize = graph.getSizeProperty("viewSrcAnchorSize")
  viewTexture = graph.getStringProperty("viewTexture")
  viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
  viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")
  viewMCLMetric = graph.getDoubleProperty("viewMCLMetric")
  viewMean = graph.getDoubleProperty("viewMean")
  viewStd = graph.getDoubleProperty("viewStd")

  """Array contenant les valeurs d'expression pour chaque temps d'expérience"""
  TP = [tp1_s, tp2_s, tp3_s, tp4_s, tp5_s, tp6_s, tp7_s, tp8_s, tp9_s, tp10_s, tp11_s, tp12_s, tp13_s, tp14_s, tp15_s, tp16_s, tp17_s]
  
  """creation d'un sous-graphe pour représenter la carte de chaleur"""
  heatmap= graph.addSubGraph("heatmap")
  
  
  """Appels des differentes fonctions"""
  Delete(graph,TP)
  Filter(graph,TP,viewMetric,Positive,Negative)
  Graphism(graph,viewShape,viewColor,viewLabel,Locus,viewBorderColor,viewBorderWidth,viewSize,Negative,Positive)
  Algorithms(graph,viewLayout,viewMCLMetric,viewMetric)
  Stats(graph, TP, viewMean, viewStd)
  Heatmap(heatmap,graph,TP, viewColor, viewSize, viewMCLMetric, viewMean, viewStd)
  
  
def Delete(graph,TP): 
  """ Permet de filtrer les gènes sans valeurs d'interactions """ 
  for n in graph.getNodes():
    if (TP[0][n] == 0 and TP[1][n] == 0 and TP[2][n] == 0) :
      graph.delNode(n)

  
    
def Filter(graph,TP,viewMetric, Positive, Negative):
  """ Permet de créer un graphe complet à partir du graphe de base
      Etablie une corrélation de Person entre chaque noeuds pour
      filtrer les arêtes peu corrélées """
  for S in graph.getNodes(): #S = Source
    for T in graph.getNodes(): # T = Target
      if S != T and S.id < T.id:
        ArrayS, ArrayT = [], []
        for i in range(1,17):
          ArrayS.append(TP[i][S])
          ArrayT.append(TP[i][T])
        P = scipy.stats.pearsonr(ArrayS,ArrayT) #result : (distance calculee, p-value)
        e = graph.existEdge(S,T,False)
        if (abs(P[1]) <= 0.01):
          if e.isValid():
            graph.delEdge(e)
        else:
          realP = 1/abs(P[0])  
          if (P[0] >= 0):     
            if e.isValid():
              viewMetric[e] = realP
              Positive[e] = True
            else:
              e=graph.addEdge(S,T)
              viewMetric[e] = realP
              Positive[e]= True
          else:     
            if e.isValid():
              viewMetric[e] = realP
              Negative[e] = True
            else:
              e=graph.addEdge(S,T)
              viewMetric[e] = realP
              Negative[e] = True
          
  
def Graphism(graph,viewShape,viewColor,viewLabel,Locus,viewBorderColor,viewBorderWidth,viewSize,Negative,Positive):
  """ Prétraitement colorant les arêtes et les noeuds """
  for n in graph.getEdges():
    viewShape[n] = tlp.EdgeShape.BezierCurve
    if (Negative[n]==False and Positive[n]==True):
      viewColor[n] = tlp.Color.Jade
    if (Negative[n]==True and Positive[n]==False):
      viewColor[n] = tlp.Color.Salmon
    if (Negative[n]==True and Positive[n]==True):
      viewColor[n] = tlp.Color.Lilac
    if (Negative[n]==False and Positive[n]==False):
      viewColor[n] = tlp.Color.Black
      
  for n in graph.getNodes():
    viewLabel[n] = Locus[n]
    viewShape[n] = tlp.NodeShape.RoundedBox
    viewColor[n] = tlp.Color.BlueGreen
    viewBorderColor[n]= tlp.Color.Black
    viewBorderWidth[n] = 3
    viewSize[n] = tlp.Size(1,0.3,1)
    if (n.id == 1307):
      viewLabel[n] = " LacZ "
      viewColor[n] = tlp.Color.Red
      viewSize[n] = tlp.Size(20,8,20)
  
def Stats(graph, TP, viewMean, viewStd):
  """ Permet de récupérer des valeurs statistiques sur les valeurs d'intéraction des gènes
      Moyenne, écart-type """
  for n in graph.getNodes():
    viewMean[n] = 0
    for i in range(len(TP)):
      viewMean[n] += TP[i][n]
    viewMean[n] /= len(TP)
    for i in range(len(TP)):
      viewStd[n] +=  (viewMean[n] - TP[i][n])**2
    viewStd[n] /= len(TP)
    viewStd[n] = math.sqrt(viewStd[n])

def Algorithms(graph,viewLayout,viewMCLMetric,viewMetric):
  """ Permet l'implémentation de l'algorithme de dessin et de cluster """
  fm3pParams = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
  graph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout, fm3pParams)
  
  #params = tlp.getDefaultPluginParameters('Fast Overlap Removal', graph)
  #graph.applyLayoutAlgorithm('Fast Overlap Removal', params)
  
  params = tlp.getDefaultPluginParameters("Edge bundling", graph)
  graph.applyAlgorithm("Edge bundling", params)

  MCLparams = tlp.getDefaultPluginParameters("MCL Clustering", graph)
  MCLparams['weights'] = viewMetric
  graph.applyDoubleAlgorithm("MCL Clustering", viewMCLMetric, MCLparams)
  
  EVparams = tlp.getDefaultPluginParameters("Equal Value", graph)
  EVparams['Property'] = viewMCLMetric
  graph.applyAlgorithm("Equal Value", EVparams)
  


def Heatmap(heatmap,graph,TP, viewColor, viewSize, viewMCLMetric, viewMean, viewStd):
  """ Permet la création de la heatmap """
  maxCluster = 0
  count=1.0
  
  for n in graph.getNodes():
    if viewMCLMetric[n] > maxCluster:
      maxCluster = int(viewMCLMetric[n])
  
  for j in range(maxCluster):
    for n in graph.getNodes():
      if viewMCLMetric[n] == j :            
        for i in range(len(TP)):
          tpValue = (TP[i][n] - viewMean[n])/viewStd[n]
          if (tpValue >= 0):
            heatmap.addNode({"viewLayout":tlp.Coord(i,count/20.0,0),  "viewColor":tlp.Color(0,int(255/2*tpValue),0), "viewSize": tlp.Size(1,0.05,1), "viewMCLMetric" : j})
          else:
            heatmap.addNode({"viewLayout":tlp.Coord(i,count/20.0,0),  "viewColor":tlp.Color(-int(255/2*tpValue),0,0), "viewSize": tlp.Size(1,0.05,1), "viewMCLMetric" : j})
        count+=1
      

