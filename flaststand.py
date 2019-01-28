import sys
import salome
salome.salome_init()
import math
import itertools

#########################################
#
# Simplified Furnace Parameters
#
#########################################

L = 16.0     # Length of the furnace
W = 4.20     # Width of the furnace
H = 1.60     # Height of the furnace
exhaustW = 3.00	    # Width of the entrance of the exhaust duct
exhaustL = 1.05	    # Length of the entrance of the exhaust duct
ro = 0.25   # Diameter of the tile exit surface
ri = 0.07   # Effective diameter of the flue gas jet
s = 1.40	# Spacing between the burners
t = 0.03 	# Size of inner box for meshing

#########################################
#
# Furnace Geometry
#
#########################################

import GEOM
from salome.geom import geomBuilder
import SALOMEDS

geompy = geomBuilder.New(salome.myStudy)
gg = salome.ImportComponentGUI("GEOM")

# Make Coordinate System
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

#------------------ Add to study
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )


#########################################
#
# Function for Furnace Surface
#
#########################################

def constructSurface(K):

    # Make Burners
    #-------------------

    #Make first burner
    centerVertex = geompy.MakeVertex(s/2.0, s/2.0, K)
    centerVector = geompy.MakeVector(centerVertex, geompy.MakeVertex(s/2, s/2, K+1))
    burnerICircle = geompy.MakeRotation(geompy.MakeCircle(centerVertex, None, ri), centerVector, math.pi/4.0)
    burnerOCircle = geompy.MakeRotation(geompy.MakeCircle(centerVertex, None, ro), centerVector, math.pi/4.0)

    # Make hexahedral divisions
    outerFace = geompy.MakeTranslationTwoPoints(geompy.MakeFaceHW(s, s, 1), O, centerVertex)
    innerFace = geompy.MakeTranslationTwoPoints(geompy.MakeFaceHW(t, t, 1), O, centerVertex)
    outerGridVertex = geompy.ExtractShapes(outerFace, geompy.ShapeType["VERTEX"], True)
    innerGridVertex = geompy.ExtractShapes(innerFace, geompy.ShapeType["VERTEX"], True)
    innerGridLines = geompy.ExtractShapes(innerFace, geompy.ShapeType["EDGE"], True)
    outerGridLines = geompy.ExtractShapes(outerFace, geompy.ShapeType["EDGE"], True)
    crossGridLines = [geompy.MakeLineTwoPnt(outerGridVertex[i], innerGridVertex[i]) for i in range(0,4)]

    compoundBurner = geompy.MakeCompound([burnerICircle, burnerOCircle]+innerGridLines+outerGridLines+crossGridLines)

    # Make rest of the burners
    burners = geompy.MakeMultiTranslation2D(compoundBurner, OX, s, 8, OY, s, 3)
    burnerLines = geompy.ExtractShapes(burners, geompy.ShapeType["EDGE"], True)

    #---------------------- Add to study
    geompy.addToStudy(centerVertex, str(K)+' center')
    geompy.addToStudy(burnerICircle, str(K)+' Inner Burner')
    geompy.addToStudy(burnerOCircle, str(K)+' Outer Burner')
    geompy.addToStudy(outerFace, str(K)+' Outer Grid Square')
    for i,j in zip(outerGridVertex, range(1,5)):
        geompy.addToStudyInFather(outerFace, i, str(K)+' OGS Vertex '+str(j))
    for i,j in zip(outerGridLines, range(1,5)):
        geompy.addToStudyInFather(outerFace, i, str(K)+' OGS Line '+str(j))
    geompy.addToStudy(innerFace, str(K)+' Inner Grid Square')
    for i,j in zip(innerGridVertex, range(1,5)):
        geompy.addToStudyInFather(innerFace, i, str(K)+' IGS Vertex '+str(j))
    for i,j in zip(innerGridLines, range(1,5)):
        geompy.addToStudyInFather(innerFace, i, str(K)+' IGS Line '+str(j))
    for i in range(0,4):
        geompy.addToStudy(crossGridLines[i], str(K)+' Cross Grid Line '+str(i+1))
    geompy.addToStudy(compoundBurner, str(K)+' Burner')
    geompy.addToStudy(burners, str(K)+' AllBurners')
    for i in burnerLines:
        geompy.addToStudyInFather(burners, i, str(K)+' Grid Line '+str(1+burnerLines.index(i)))


    # Make exhaust and wall
    #-----------------------

    # Make exhaust
    exhaustVertex = [geompy.MakeVertex(L, (W+exhaustW)/2.0, K),
        geompy.MakeVertex(L, (W-exhaustW)/2.0, K),
        geompy.MakeVertex(L-exhaustL, (W-exhaustW)/2.0, K),
        geompy.MakeVertex(L-exhaustL, (W+exhaustW)/2.0, K)]
    exhaustLine = [geompy.MakeLineTwoPnt(exhaustVertex[1], exhaustVertex[2]),
        geompy.MakeLineTwoPnt(exhaustVertex[3], exhaustVertex[0])]

    # Make walls
    wallVertex = [geompy.MakeVertex(L-exhaustL, W, K), geompy.MakeVertex(L-exhaustL, 0, K),
        geompy.MakeVertex(s*8, s, K), geompy.MakeVertex(s*8, s*2, K)]
    wallLine = [geompy.MakeLineTwoPnt(wallVertex[0], wallVertex[1]),
        geompy.MakeLineTwoPnt(wallVertex[2], exhaustVertex[2]),
        geompy.MakeLineTwoPnt(wallVertex[3], exhaustVertex[3])]

    # ---------------------------- Add to study
    geompy.addToStudy( exhaustLine[0], 'Exhaust Line 1' )
    geompy.addToStudyInFather( exhaustLine[0], exhaustVertex[1], 'Exhaust Vertex 2' )
    geompy.addToStudyInFather( exhaustLine[0], exhaustVertex[1], 'Exhaust Vertex 3' )
    geompy.addToStudy( exhaustLine[1], 'Exhaust Line 2' )
    geompy.addToStudyInFather( exhaustLine[1], exhaustVertex[3], 'Exhaust Vertex 4' )
    geompy.addToStudyInFather( exhaustLine[1], exhaustVertex[0], 'Exhaust Vertex 1' )
    geompy.addToStudy(wallLine[0], 'Wall Line 1')
    geompy.addToStudyInFather( wallLine[0], wallVertex[0], 'Wall Vertex 1' )
    geompy.addToStudyInFather( wallLine[0], wallVertex[1], 'Wall Vertex 2' )
    geompy.addToStudy(wallLine[1], 'Wall Line 2')
    geompy.addToStudyInFather( wallLine[1], wallVertex[2], 'Wall Vertex 3' )
    geompy.addToStudy(wallLine[2], 'Wall Line 3')
    geompy.addToStudyInFather( wallLine[2], wallVertex[3], 'Wall Vertex 4' )


    # Make the surface
    centerFurnace = geompy.MakeVertex(L/2.0, W/2.0, K)
    baseSurface = geompy.MakeTranslationTwoPoints(geompy.MakeFaceHW(L, W, 1), O, centerFurnace)

    # Partition the upper face into burners and exhaust
    Partition = geompy.MakePartition([baseSurface], burnerLines+wallLine+exhaustLine, [], [], geompy.ShapeType["FACE"], 0, [], 0)
    surface = geompy.ExtractShapes(Partition, geompy.ShapeType["FACE"], True)

    # ---------------------------- Add to study
    geompy.addToStudy( centerFurnace, str(K)+' Furnace center' )
    geompy.addToStudy( baseSurface, str(K)+' face' )
    geompy.addToStudy( Partition, str(K)+' face Partition' )
    for i in range(len(surface)):
        geompy.addToStudyInFather( Partition, surface[i], str(K)+' surface '+str(i+1) )

    return surface




#---------------------------
# Make Furnace
#---------------------------

# Make surface
topSurface = constructSurface(H)
bottomSurface = constructSurface(0)

# Connect with hexahedrons
furnaceSolid = [geompy.MakeHexa2Faces(topSurface[i], bottomSurface[i]) for i in range(len(topSurface))]
hexaSideFaces = [geompy.ExtractShapes(i, geompy.ShapeType["FACE"], True) for i in furnaceSolid]
compoundFurnace = geompy.MakeCompound(furnaceSolid)
finalFurnace = geompy.MakeGlueFaces(compoundFurnace, 1.0e-7)

# ---------------------------- Add to study
geompy.addToStudy(compoundFurnace, 'Compound Furnace')
for i in range(len(furnaceSolid)):
    geompy.addToStudyInFather(compoundFurnace, furnaceSolid[i], 'Furnace Hexa '+str(i))
    for j in range(6):
        geompy.addToStudyInFather(furnaceSolid[i], hexaSideFaces[i][j], 'Furnace Hexa '+str(i)+' Face '+str(j))
geompy.addToStudy(finalFurnace, 'Final Furnace')







#########################################
#
# Furnace Mesh
#
#########################################

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(salome.myStudy)

hexMesh = smesh.Mesh(finalFurnace)
Regular_1D = hexMesh.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(20)
Quadrangle_2D = hexMesh.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = hexMesh.Hexahedron(algo=smeshBuilder.Hexa)
isDone = hexMesh.Compute()

#----------------------- Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(hexMesh.GetMesh(), 'Hexa Mesh')


# Group Exhuast
#---------------
aCriteria = [smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[316])]
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(hexMesh.GetMesh())
exhaust = hexMesh.GroupOnFilter( SMESH.FACE, 'exhaust', aFilter )

# --------- Set name of Mesh object
smesh.SetName(exhaust, 'exhaust')


# Group Burner Tiles
#--------------------
tileId = [	[3+39*i, 7+39*i, 11+39*i, 22+39*i] for i in range(8)	] + [
          	[6+39*i, 17+39*i, 21+39*i, 32+39*i] for i in range(8)	] + [
          	[16+39*i, 27+39*i, 31+39*i, 35+39*i] for i in range(8)	]
for i in range(24):
    aCriteria = [
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[tileId[i][0]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[tileId[i][1]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[tileId[i][2]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[tileId[i][3]])	]
    aFilter = smesh.GetFilterFromCriteria(aCriteria)
    aFilter.SetMesh(hexMesh.GetMesh())
    burnerTile = hexMesh.GroupOnFilter( SMESH.FACE, 'burnerTile_'+str(i+1), aFilter )
    burnerTile.SetColor( SALOMEDS.Color( 1, 0, 0 ))
    # --------- Set name of Mesh object
    smesh.SetName(burnerTile, 'burnerTile_'+str(i+1))


# Group Burners type A
#----------------
burnerAId = [	[5+39*i, 8+39*i, 9+39*i, 10+39*i, 12+39*i] for i in range(3)	] + [
          		[15+39*i, 18+39*i, 19+39*i, 20+39*i, 23+39*i] for i in range(3)	] + [
          		[26+39*i, 28+39*i, 29+39*i, 30+39*i, 33+39*i] for i in range(3)	]
for i in range(9):
    aCriteria = [
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerAId[i][0]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerAId[i][1]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerAId[i][2]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerAId[i][3]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerAId[i][4]])]
    aFilter = smesh.GetFilterFromCriteria(aCriteria)
    aFilter.SetMesh(hexMesh.GetMesh())
    burnerA = hexMesh.GroupOnFilter( SMESH.FACE, 'burnerA_'+str(i+1), aFilter )
    burnerA.SetColor( SALOMEDS.Color( 0.4, 1, 0 ))
    # --------- Set name of Mesh object
    smesh.SetName(burnerA, 'burnerA_'+str(i+1))


# Group Burners type B
#----------------
burnerBId = [	[5+39*i, 8+39*i, 9+39*i, 10+39*i, 12+39*i] for i in range(3,8)	] + [
          		[15+39*i, 18+39*i, 19+39*i, 20+39*i, 23+39*i] for i in range(3,8)	] + [
          		[26+39*i, 28+39*i, 29+39*i, 30+39*i, 33+39*i] for i in range(3,8)	]
for i in range(15):
    aCriteria = [
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerBId[i][0]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerBId[i][1]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerBId[i][2]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerBId[i][3]],SMESH.FT_Undefined,SMESH.FT_LogicalOR),
    smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topSurface[burnerBId[i][4]])]
    aFilter = smesh.GetFilterFromCriteria(aCriteria)
    aFilter.SetMesh(hexMesh.GetMesh())
    burnerB = hexMesh.GroupOnFilter( SMESH.FACE, 'burnerB_'+str(i+1), aFilter )
    burnerB.SetColor( SALOMEDS.Color( 0, 1, 0.4 ))
    # --------- Set name of Mesh object
    smesh.SetName(burnerB, 'burnerB_'+str(i+1))


# Group Top Wall
#-----------------
setOfBurners = set([topSurface[i] for j in burnerAId for i in j]+[topSurface[i] for j in burnerBId for i in j]+[topSurface[i] for j in tileId for i in j]+[topSurface[316]])
topWallSurface = [x for x in topSurface if x not in setOfBurners]
aCriteria = []
for i in topWallSurface[:-1]:
    aCriteria.append(smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,i,SMESH.FT_Undefined,SMESH.FT_LogicalOR))
aCriteria.append(smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,topWallSurface[-1]))
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(hexMesh.GetMesh())
topWall = hexMesh.GroupOnFilter( SMESH.FACE, 'topWall', aFilter )
topWall.SetColor( SALOMEDS.Color( 0, 0, 0.3 ))

# --------- Set name of Mesh object
smesh.SetName(topWall, 'topWall')


# Group Bottom Wall
#-------------------
aCriteria = []
for i in bottomSurface[:-1]:
    aCriteria.append(smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,i,SMESH.FT_Undefined,SMESH.FT_LogicalOR))
aCriteria.append(smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,bottomSurface[-1]))
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(hexMesh.GetMesh())
bottomWall = hexMesh.GroupOnFilter( SMESH.FACE, 'bottomWall', aFilter )
bottomWall.SetColor( SALOMEDS.Color( 0, 0, 0.9 ))

# --------- Set name of Mesh object
smesh.SetName(bottomWall, 'bottomWall')


# Group Side Wall
#-------------------
wallId = [[4+39*i,1] for i in range(8)]+[[312,3],[315,1],[315,5],[316,5],[317,5],[317,4],[313,4],[2,0],[1,0],[0,0]]+[[34+39*i,4] for i in range(8)]
aCriteria = []
for i in wallId[:-1]:
	aCriteria.append(smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,hexaSideFaces[i[0]][i[1]],SMESH.FT_Undefined,SMESH.FT_LogicalOR))
aCriteria.append(smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGeom,SMESH.FT_Undefined,hexaSideFaces[wallId[-1][0]][wallId[-1][1]]))
aFilter = smesh.GetFilterFromCriteria(aCriteria)
aFilter.SetMesh(hexMesh.GetMesh())
sideWalls = hexMesh.GroupOnFilter( SMESH.FACE, 'sideWalls', aFilter )
sideWalls.SetColor( SALOMEDS.Color( 0, 0, 0.6 ) )

# --------- Set name of Mesh object
smesh.SetName(sideWalls, 'sideWalls')



if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
