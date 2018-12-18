# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/lintj/Desktop')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

box = geompy.MakeBoxDXDYDZ(16, 5.1, 1.5)
geompy.addToStudy( box, 'box' )
boxFaces = geompy.ExtractShapes(box, geompy.ShapeType["FACE"], True)
for i in range(6):
    geompy.addToStudyInFather( box, boxFaces[i], 'boxface_'+str(i))

iCircle = [] 
oCircle = []

for i in range(0,8):
    Vertex_1 = geompy.MakeVertex(0.48+1.40*i, 1.15, 1.5)
    iCircle.append(geompy.MakeCircle(Vertex_1, None, 0.0865))
    oCircle.append(geompy.MakeCircle(Vertex_1, None, 0.2875))

    Vertex_2 = geompy.MakeVertex(0.48+1.40*i, 2.55, 1.5)
    iCircle.append(geompy.MakeCircle(Vertex_2, None, 0.0865))
    oCircle.append(geompy.MakeCircle(Vertex_2, None, 0.2875))

    Vertex_3 = geompy.MakeVertex(0.48+1.40*i, 3.95, 1.5)
    iCircle.append(geompy.MakeCircle(Vertex_3, None, 0.0865))
    oCircle.append(geompy.MakeCircle(Vertex_3, None, 0.2875))

Vertex_39 = geompy.MakeVertex(16, 4.125, 1.5)
Vertex_40 = geompy.MakeVertex(16, 0.975, 1.5)
Vertex_41 = geompy.MakeVertex(15.05, 0.975, 1.5)
Vertex_42 = geompy.MakeVertex(15.05, 4.125, 1.5)
Line_1 = geompy.MakeLineTwoPnt(Vertex_42, Vertex_41)
Line_2 = geompy.MakeLineTwoPnt(Vertex_41, Vertex_40)
Line_3 = geompy.MakeLineTwoPnt(Vertex_39, Vertex_42)

for i in range(24):
    geompy.addToStudyInFather( box, iCircle[i], 'iCircle_'+str(i+1) )
    geompy.addToStudyInFather( box, oCircle[i], 'oCircle_'+str(i+1) )
geompy.addToStudyInFather( box, Line_1, 'Line_1' )
geompy.addToStudyInFather( box, Line_2, 'Line_2' )
geompy.addToStudyInFather( box, Line_3, 'Line_3' )

Partition_1 = geompy.MakePartition([boxFaces[3]], iCircle+oCircle+[Line_1, Line_2, Line_3], [], [], geompy.ShapeType["FACE"], 0, [], 0)
topFurnace = geompy.ExtractShapes(Partition_1, geompy.ShapeType["FACE"], True)
geompy.addToStudy( Partition_1, 'Partition_1' )
for i in range(len(topFurnace)):
    geompy.addToStudyInFather( Partition_1, topFurnace[i], 'topFace_'+str(i+1) )

furnace = geompy.MakeSolidFromConnectedFaces([boxFaces[0],boxFaces[1],boxFaces[2],boxFaces[4],boxFaces[5]]+topFurnace,True)
geompy.addToStudy( furnace, 'furnace' )

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
