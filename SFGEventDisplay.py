from ROOT import gSystem
from ROOT import gStyle
from ROOT import TFile
from ROOT import TArrayI
from ROOT import TArrayD
from ROOT import TTree
from ROOT import std
from ROOT import TMath
from ROOT import TVector3
from ROOT import TH1D
from ROOT import TH2D
from ROOT import TCanvas
import numpy as np
import sys
import math
from array import array
import re
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def SetPlotSettings(plot,option):

  plot.GetXaxis().SetLabelSize(0.05)
  plot.GetXaxis().SetTitleSize(0.05)
  plot.GetYaxis().SetLabelSize(0.05)
  plot.GetYaxis().SetTitleSize(0.05)
  plot.GetYaxis().SetTitleOffset(1.3)

  plot.GetZaxis().SetLabelSize(0.04)
  
  #plot.GetXaxis().SetTitle("Z axis position / mm")
  if(option == 1):
    plot.GetXaxis().SetTitle("Z axis position")
    plot.GetYaxis().SetTitle("Y axis position")
    plot.SetTitle("YZ view")
  elif(option == 2):
    plot.GetXaxis().SetTitle("Z axis position")
    plot.GetYaxis().SetTitle("X axis position")
    plot.SetTitle("XZ view")
  elif(option == 3):
    plot.GetXaxis().SetTitle("X axis position")
    plot.GetYaxis().SetTitle("Y axis position")
    plot.SetTitle("XY view")

# --------------------------------
# Below are the main codes
# --------------------------------

# Try to load the project shared library. Might already be built.
buildProject = False
if (gSystem.Load("SFGAnalysisProject/SFGAnalysisProject.so") < 0):
        buildProject = True
        
# --------------------------------
# Get the ROOT TFile we will read
# --------------------------------
file = TFile(sys.argv[1],"OLD")

print ("Input file launched:" + sys.argv[1])

# Build the project shared library if it's needed. This loads it too.
if (buildProject):
        file.MakeProject("SFGAnalysisProject","ND::TSFGReconModule","recreate++")

# Get the SFG reconstruction tree. There are other trees too.
sfgTree = file.Get("ReconDir/SFG")
trajTree = file.Get("TruthDir/Trajectories")
verTree = file.Get("TruthDir/Vertices")

nEntries = sfgTree.GetEntries()

print ("Entries = ")
print (nEntries)

# --------------------------------
# Create output file, tree and variables
# --------------------------------
out_file = TFile(sys.argv[2],"recreate")

print ("Output file settled: " + sys.argv[2])

# Define some variables there

gStyle.SetOptStat(0)

n_itr = 0

file.cd()

# --------------------------------
# Loop over all events in the tree
# --------------------------------

# for entry1 in trajTree:
#   trajTree.GetEntry(n_itr)
#   n_itr+=1
#   print ("EventID = " )
#   print (trajTree.EventID)
#   fig1 = plt.figure()
#   ax = plt.axes(projection='3d')
#   for trj in trajTree.Trajectories:
#     print ("PDG = " )
#     print ( trj.PDG)
#     print ("ID = " )
#     print ( trj.ID)
#     tmp = trj.InitPosition
#     n_hit = 0;
#     hit_posx = []
#     hit_posy = []
#     hit_posz = []
#     hit_posx.append(tmp.X())
#     hit_posy.append(tmp.Y())
#     hit_posz.append(tmp.Z())
#     for pt in trj.Points:
#       hit_posx.append(pt.PositionX)
#       hit_posy.append(pt.PositionY)
#       hit_posz.append(pt.PositionZ)
#     tmp = trj.FinalPosition
#     hit_posx.append(tmp.X())
#     hit_posy.append(tmp.Y())
#     hit_posz.append(tmp.Z())
#     #ax.scatter(hit_posz,hit_posx,hit_posy,linestyle='-',linewidth=4,marker='s',edgecolors='b',label='Hit')
#     ax.plot(hit_posz,hit_posx,hit_posy,linestyle='-',linewidth=4,marker='s',label='{}'.format(trj.PDG))
#     ax.legend(fontsize=10,numpoints=1)
#     ax.set_xlabel('Z axis position')
#     ax.set_ylabel('X axis position')
#     ax.set_zlabel('Y axis position')
#     #ax.set_xlim3d(0,184)
#     #ax.set_ylim3d(0,192)
#     #ax.set_zlim3d(0,56)  
#     #plt.savefig('./SFGevent3D/event_%i_hit' % (index))
#   plt.show()

for entry in sfgTree:

  index = entry.EventID
 
  # if(index > 567800000):
  #   break
  
#   if(index != 25):
#     print ("Not 25!")
#     continue
  
  if(entry.NAlgoResults == 0 or entry.NParticles == 0):
    print ("No results!")
    continue

  # Create the 2D histograms to draw event YZ, XZ, XY views
  event_yz = TH2D("event_yz","event_yz",184,0,184,56,0,56)
  SetPlotSettings(event_yz,1)

  event_xz = TH2D("event_xz","event_xz",184,0,184,192,0,192)
  SetPlotSettings(event_xz,2)
  
  event_xy = TH2D("event_xy","event_xy",192,0,192,56,0,56)
  SetPlotSettings(event_xy,3)

  # Loop over all hits
  n_hit = 0;
  hit_posx = []
  hit_posy = []
  hit_posz = []
  
  for hit in entry.AlgoResults[0].Hits:
  #for hit in entry.TrueHits:
  
    temp = '{:032b}'.format(entry.Hits[hit].GeomId)
    #temp = '{:032b}'.format(entry.TrueHits[n_hit].GeomId)
    hitx = int(temp[10:18],2)
    hity = int(temp[18:24],2)
    hitz = int(temp[24:32],2)
    
    hitcharge = entry.Hits[hit].Charge
    #hitcharge = entry.TrueHits[n_hit].Charge
    
    hit_posx.append(hitx)
    hit_posy.append(hity)
    hit_posz.append(hitz)
    
    event_yz.SetBinContent(hitz,hity,hitcharge)
    event_xz.SetBinContent(hitz,hitx,hitcharge)
    event_xy.SetBinContent(hitx,hity,hitcharge)
    
    n_hit += 1
 
  # Loop over all nodes
  n_node = 0
  node_posx = []
  node_posy = []
  node_posz = []

  # for par in entry.AlgoResults[0].Particles:
  
  #   trk_count = 0
  
  #   # node_posx = []
  #   # node_posy = []
  #   # node_posz = []

  #   for trk in entry.Particles[par].Tracks:
 
  #     if(trk_count != 0):
  #       continue
      
  #     trk_count += 1
      
  #     for node in entry.Tracks[trk].Nodes:
      
  #       temp = entry.Nodes[node].Position
  #       node_posx.append(temp.X())
  #       node_posy.append(temp.Y())
  #       node_posz.append(temp.Z())
 
  # Draw plots
    
  # 3D plot
  fig1 = plt.figure()
  ax = plt.axes(projection='3d')
  ax.scatter(hit_posz,hit_posx,hit_posy,marker='s',edgecolors='b',label='Hit')
  ax.legend(fontsize=10,numpoints=1)
  ax.set_xlabel('Z axis position')
  ax.set_ylabel('X axis position')
  ax.set_zlabel('Y axis position')
  #ax.set_xlim3d(0,184)
  #ax.set_ylim3d(0,192)
  #ax.set_zlim3d(0,56)  
  fig1.savefig('./SFGevent3D/event_%i_hit.pdf' % (index))
  # plt.show(block=False)
  # plt.pause(0)
  plt.close(fig1)

  colors=['b','g','r','c','m','y','k']
  fig2 = plt.figure()
  ax = plt.axes(projection='3d')
  unused_posx = []
  unused_posy = []
  unused_posz = []
  for hit in entry.Hits:
    temp = '{:032b}'.format(hit.GeomId)
    hitx = (int(temp[10:18],2)+0.5)*10.27-985.92
    hity = (int(temp[18:24],2)+0.5)*10.27+30-287.56
    hitz = (int(temp[24:32],2)+0.5)*10.27-2888.78
    unused_posx.append(hitx)
    unused_posy.append(hity)
    unused_posz.append(hitz)
  if len(unused_posx)>0:
    ax.scatter(unused_posz,unused_posx,unused_posy,marker='o',color="None",edgecolors='grey',label='Unused',s=1)
    ax.plot([unused_posz[0]],[unused_posx[0]],[unused_posy[0]],marker='o',linestyle="none",color='grey',label='Hits')
  for par in entry.AlgoResults[0].Particles:
  
    trk_count = 0
  
    node_posx = []
    node_posy = []
    node_posz = []

    for trk in entry.Particles[par].Tracks:
 
      if(trk_count != 0):
        continue
      
      trk_count += 1
      
      for node in entry.Tracks[trk].Nodes:
      
        temp = entry.Nodes[node].Position
        node_posx.append(temp.X())
        node_posy.append(temp.Y())
        node_posz.append(temp.Z())

    ax.scatter(node_posz,node_posx,node_posy,marker='o',color=colors[n_node%7],label='Node{}'.format(n_node))
    ax.plot([node_posz[0]],[node_posx[0]],[node_posy[0]],marker='o',linestyle="none",color=colors[n_node%7],label='Node{}'.format(n_node))
    n_node+=1
  #ax.scatter(node_posz,node_posx,node_posy,marker='o',edgecolors='b',label='Node')
  #legend1 = ax.legend(fontsize=10,numpoints=1,loc='lower right')
  #plt.gca().add_artist(legend1)
  #ax.set_xlim3d(0,184)
  #ax.set_ylim3d(0,192)
  #ax.set_zlim3d(0,56)  
  #if entry.AlgoResults[0].UnusedHitCluster > 0:
  xmin, xmax = ax.get_xlim3d()
  ymin, ymax = ax.get_ylim3d()
  zmin, zmax = ax.get_zlim3d()
  #plt.savefig('./SFGevent3D/event_%i_node' % (index))
  pdgId=[211,-211,11,-11,13,-13,22,2212,2112,999]
  pdgName=['pi+','pi-','e-','e+','mu-','mu+','gamma','p','n','other']
  pdgColor=['blue','cyan','green','lime','red','darkred','yellow','black','grey','violet']
  pdgUse=[0]*len(pdgId)
  pdgPlot=[]
  pdgPlotLabel=[]
  trajTree.GetEntry(n_itr)
  for trj in trajTree.Trajectories:
    # if (trj.PDG==2212 or trj.PDG==2112):
    #   continue
    tmp = trj.InitPosition
    hit_posx = []
    hit_posy = []
    hit_posz = []
    hit_posx.append(tmp.X())
    hit_posy.append(tmp.Y())
    hit_posz.append(tmp.Z())
    for pt in trj.Points:
      hit_posx.append(pt.PositionX)
      hit_posy.append(pt.PositionY)
      hit_posz.append(pt.PositionZ)
    tmp = trj.FinalPosition
    hit_posx.append(tmp.X())
    hit_posy.append(tmp.Y())
    hit_posz.append(tmp.Z())
    idx=-1
    if trj.PDG in pdgId:
      idx=pdgId.index(trj.PDG)
    #tmp = ax.plot(hit_posz,hit_posx,hit_posy,linestyle='-',marker='s',color=pdgColor[idx],label='{}'.format(pdgName[idx]))
    if pdgUse[idx]==0:
      ax.plot(hit_posz,hit_posx,hit_posy,linestyle='-',marker='s',color=pdgColor[idx],label='{}'.format(pdgName[idx]))
      # pdgPlot.append(tmp)
      # pdgPlotLabel.append(pdgName[idx])
      pdgUse[idx]=1
    else:
      ax.plot(hit_posz,hit_posx,hit_posy,linestyle='-',marker='s',color=pdgColor[idx],label='_nolegend_')
    #   ax.plot(hit_posz,hit_posx,hit_posy,linestyle='-',marker='s',color=pdgColor[-1],label='{}'.format(pdgName[-1]))
  legend2 = ax.legend(fontsize=10,numpoints=1,loc='upper left')
  # ax.add_artist(legend1)
  # ax.add_artist(legend2)
  ax.set_xlabel('Z axis position')
  ax.set_ylabel('X axis position')
  ax.set_zlabel('Y axis position')
  # plt.ylim(-1000, 1000)
  # plt.xlim(-3000,-1000)
  # ax.set_xlim3d(-3000,-1000)
  # ax.set_ylim3d(-1000, 1000)
  # ax.set_zlim3d(-300,300)
  ax.set_xlim3d(xmin,xmax)
  ax.set_ylim3d(ymin,ymax)
  ax.set_zlim3d(zmin,zmax)
  fig2.savefig('./SFGevent3D/event_%i_truth.pdf' % (index))
  # plt.show(block=False)
  # plt.pause(0)
  plt.close(fig2)
  
  # # 2D plot
  # plot_yz = TCanvas("plot_yz","plot_yz",800,600)
  # plot_yz.SetRightMargin(0.15)
  # plot_yz.SetLeftMargin(0.15)
  # plot_yz.SetBottomMargin(0.13)
  # plot_yz.cd()
  # event_yz.SetContour(99)
  # event_yz.Draw("colz")
  # plot_yz.Update()

  # plot_xz = TCanvas("plot_xz","plot_xz",800,600)
  # plot_xz.SetRightMargin(0.15)
  # plot_xz.SetLeftMargin(0.15)
  # plot_xz.SetBottomMargin(0.13)
  # plot_xz.cd()
  # event_xz.SetContour(99)
  # event_xz.Draw("colz")
  # plot_xz.Update()
  
  # plot_xy = TCanvas("plot_xy","plot_xy",800,600)
  # plot_xy.SetRightMargin(0.15)
  # plot_xy.SetLeftMargin(0.15)
  # plot_xy.SetBottomMargin(0.13)
  # plot_xy.cd()
  # event_xy.SetContour(99)
  # event_xy.Draw("colz")
  # plot_xy.Update()   
 
  # # Save event plot in output file
  # out_file.cd()
  # plot_yz.Write("event_yz_%i" % (index))
  # plot_xz.Write("event_xz_%i" % (index))
  # plot_xy.Write("event_xy_%i" % (index))
  # file.cd()
  
  n_itr += 1
    
# --------------------------------
# Finish
# -------------------------------- 

out_file.Close()

file.Close()

print ("ALL FINISHED!")
