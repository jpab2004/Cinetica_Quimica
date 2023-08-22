#===============================================================================================#
#                                             Setup                                             #
#===============================================================================================#
# Libraries
from itertools import combinations
from random import uniform
from pickle import load
from math import isnan
from vpython import *

# Preparing Scene
scene.delete()
sceneBuffer = 1
sceneWidth = 900 * sceneBuffer
sceneHeight = 900 * sceneBuffer
scene = canvas(width=sceneHeight, height=sceneHeight, align='left')
scene.background = vector(0, 0, 0)
scene.append_to_caption("<div id='fps'/>")

# Elements
file = 'Assets/elementsHashTable.pickle'
with open(file, 'rb') as f:
    elements = load(f)



#===============================================================================================#
#                                        Common Functions                                       #
#===============================================================================================#
def RGB2VEC(r, g, b): 
    return vector(r/255, g/255, b/255)

def createWalls():
    if solidWalls:
        s2 = 2 * side - thickness
        s3 = 2 * side + thickness

        wallE = box(pos=vector(side, 0, 0), size=vector(thickness, s2, s2), color=vector(1, 0, 0))
        wallW = box(pos=vector(-side, 0, 0), size=vector(thickness, s2, s2), color=vector(1, 0, 0))

        wallN = box(pos=vector(0, side, 0), size=vector(s3, thickness, s2), color=vector(0, 1, 0))
        wallS = box(pos=vector(0, -side, 0), size=vector(s3, thickness, s2), color=vector(0, 1, 0))
        
        wallC = box(pos=vector(0, 0, -side), size=vector(s3, s3, thickness), color=vector(0, 0, 1))
    else:
        edgeColor = vector(1, 0, 1)
        edgeRadius = .2
        edgeLength = 2*side

        if overallPretty:
            vertex000 = sphere(pos=vector(-side, -side, -side), radius=edgeRadius, color=edgeColor)
            vertex001 = sphere(pos=vector(-side, -side, side), radius=edgeRadius, color=edgeColor)
            vertex010 = sphere(pos=vector(-side, side, -side), radius=edgeRadius, color=edgeColor)
            vertex011 = sphere(pos=vector(-side, side, side), radius=edgeRadius, color=edgeColor)
            vertex100 = sphere(pos=vector(side, -side, -side), radius=edgeRadius, color=edgeColor)
            vertex101 = sphere(pos=vector(side, -side, side), radius=edgeRadius, color=edgeColor)
            vertex110 = sphere(pos=vector(side, side, -side), radius=edgeRadius, color=edgeColor)
            vertex111 = sphere(pos=vector(side, side, side), radius=edgeRadius, color=edgeColor)

        edgeSB = cylinder(pos=vector(-side, -side, -side), axis=vector(1, 0, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeNB = cylinder(pos=vector(-side, side, -side), axis=vector(1, 0, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeSF = cylinder(pos=vector(-side, -side, side), axis=vector(1, 0, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeNF = cylinder(pos=vector(-side, side, side), axis=vector(1, 0, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)

        edgeSW = cylinder(pos=vector(-side, -side, -side), axis=vector(0, 0, 1), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeNW = cylinder(pos=vector(-side, side, -side), axis=vector(0, 0, 1), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeSE = cylinder(pos=vector(side, -side, -side), axis=vector(0, 0, 1), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeNE = cylinder(pos=vector(side, side, -side), axis=vector(0, 0, 1), length=edgeLength, color=edgeColor, radius=edgeRadius)

        edgeWB = cylinder(pos=vector(-side, -side, -side), axis=vector(0, 1, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeWF = cylinder(pos=vector(-side, -side, side), axis=vector(0, 1, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeEB = cylinder(pos=vector(side, -side, -side), axis=vector(0, 1, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)
        edgeEF = cylinder(pos=vector(side, -side, side), axis=vector(0, 1, 0), length=edgeLength, color=edgeColor, radius=edgeRadius)

    return

def generateTheoryCurve(e, nParticlesTheory):
    theory = gcurve(color=elements[e]['color'], label=elements[e]['name-pt'])
    mass = generateMass(e)


    deltaV = 10
    for v in range(0, maxVel+deltaV, deltaV):
        alpha = (dv**2/deltaV) * nParticlesTheory
        first = (mass / (2*pi*k*temperature))**1.5
        second = exp(-.5*mass*v**2 / (k*temperature))*v**2

        value = alpha * first * second
        theory.plot(v, value)

def generatePosition(d3=True):
    if randomPosition:
        if d3: return positionBuffer*vector.random()
        else: return vector(positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1), 0)
    else: return vector(0, 0, 0)

def generateVelocity(particleMass, d3=True):
    averageKinecticMomentum = sqrt(2*particleMass*1.5*k*temperature)

    if d3:
        psi = radians(uniform(0, 360))
        theta = radians(uniform(0, 360))

        x = averageKinecticMomentum * cos(theta) * sin(psi)
        y = averageKinecticMomentum * sin(theta) * sin(psi)
        z = averageKinecticMomentum * cos(psi)

        return vector(x, y, z)
    else:
        theta = radians(uniform(0, 360))

        x = averageKinecticMomentum * cos(theta)
        y = averageKinecticMomentum * sin(theta)

        return vector(x, y, 0)
    
def generateMass(e):
    return elements[e]['mass']*1E-3/6e23

def generateParticle(e, d3=True):
    particleRadius = getRadii(e)
    particleMass = generateMass(e)
    particleColor = elements[e]['color']
    particleVelocity = generateVelocity(particleMass, d3)/particleMass

    still = True
    particlePosition = generatePosition(d3)
    positions = list(map(getPositionAndRadius, particles))
    while (still):
        still = False
        particlePosition = generatePosition(d3)
        for pos, rad in positions:
            if (mag(pos - particlePosition)) <= (rad + particleRadius):
                still = True
                break

    particle = sphere(pos=particlePosition, radius=particleRadius, color=particleColor, make_trail=makeTrails, retain=10)
    particle.v = particleVelocity
    particle.m = particleMass
    return particle

def getRadii(e):
    f = lambda x: (exp(x/500) - .7) * radiiBuff

    if empiricalRadii:
        radii = elements[e]['radii-empirical']
    else:
        radii = elements[e]['radii-calculated']

    if isnan(radii): return f(150)    

    return f(radii)

def getMaxRadii():
    maxRadii = float('-inf')

    if empiricalRadii: radiiChosen = 'radii-empirical'
    else: radiiChosen = 'radii-calculated'

    for e in elementsToSimulate:
        if elements[e][radiiChosen] > maxRadii: maxRadii = elements[e][radiiChosen]
    
    return maxRadii

def getPositionAndRadius(p): return (p.pos, p.radius)

def startSimulation():
    global globalStart
    globalStart = True
    startSimulationButton.disabled = True

    return

def newVelocity(p1, p2):
    v1, v2 = p1.v, p2.v
    m1, m2 = p1.m, p2.m
    r1, r2 = p1.pos, p2.pos

    v1_ = v1 - ((2*m2) / (m1 + m2)) * (dot(v1 - v2, r1 - r2) / mag2(r1 - r2)) * (r1 - r2)
    v2_ = v2 - ((2*m1) / (m2 + m1)) * (dot(v2 - v1, r2 - r1) / mag2(r2 - r1)) * (r2 - r1)

    return v1_, v2_

def collision(iterator):
    for p1, p2 in combinations(iterator, 2):
        if (mag(p1.pos - p2.pos) <= p1.radius + p2.radius):
            p1.v, p2.v = newVelocity(p1, p2)
    
    return



#===============================================================================================#
#                                          3D Functions                                         #
#===============================================================================================#
def step3D(iterator):
    for p in iterator:
        p.pos += p.v*dt
        wallCollision = side - .5*thickness - p.radius

        if not (wallCollision > p.pos.x > -wallCollision):
            p.v.x *= -1
            p.pos += p.v*dt
        if not (wallCollision > p.pos.y > -wallCollision):
            p.v.y *= -1
            p.pos += p.v*dt
        if not (wallCollision > p.pos.z > -wallCollision):
            p.v.z *= -1
            p.pos += p.v*dt

    return



#===============================================================================================#
#                                          2D Functions                                         #
#===============================================================================================#
def step2D(iterator):
    for p in iterator:
        p.pos += p.v*dt
        wallCollision = side - .5*thickness - p.radius

        if not (wallCollision > p.pos.x > -wallCollision):
            p.v.x *= -1
            p.pos += p.v*dt
        if not (wallCollision > p.pos.y > -wallCollision):
            p.v.y *= -1
            p.pos += p.v*dt
    
    return



#===============================================================================================#
#                                       Running Function                                        #
#===============================================================================================#
def run(d3=True):
    if d3: runFunction = step3D
    else: runFunction = step2D

    for element, count in zip(elementsToSimulate, elementsCount):
        generateTheoryCurve(element, count)
        for _ in range(count):
            particle = generateParticle(element, d3)
            particles.append(particle)

    while(not globalStart):
        rate(15)
        continue

    i = 1

    while(globalStart):
        rate(fps)
        runFunction(particles)
        collision(particles)

        if i >= loopVerboseCount: drawHist(particles); i = 1
        i += 1
    
    return



#===============================================================================================#
#                                          Controllers                                          #
#===============================================================================================#
startSimulationButton = button(pos=scene.caption_anchor, text='Start simulation', bind=startSimulation)



#===============================================================================================#
#                                        Velocities Graph                                       #
#===============================================================================================#
def getVelocity(p): return mag(p.v)
def getHist(v):
    return int(v/dv)

def drawHist(particles):
    global bars

    vels = list(map(getVelocity, particles))
    histData = {}
    for i in list(map(getHist, vels)):
        try:
            histData[i] += 1
        except:
            histData[i] = 0

    histData = [[v*dv + .5*dv, histData[v]] if v in histData else [v*dv + .5*dv, 0] for v in range(max(histData))]
    bars.data = histData


#===============================================================================================#
#                                           Simulation                                          #
#===============================================================================================#
# Wall variables
side = 10
thickness = .5

# Prettier
overallPretty = False
makeTrails = False
solidWalls = False

# Particle variables
particles = []
radiiBuff = .2
positionBuffer = 8
randomPosition = True
empiricalRadii = True

# Elements
elementsToSimulate = [2]
elementsCount = [600]
nParticles = sum(elementsCount)

# Velocity graph variables
dv = 100
maxVel = 6000
graphWidth = 800
histData = [(v*dv + .5*dv, 0) for v in range(int(maxVel/dv))]
velGraph = graph(title='Velocidade das partículas na simulação', xtitle='Velocidade (m/s)', xmax=maxVel, ymax=nParticles/4,
                 ytitle='Número de Partículas', fast=False, width=800, align='left', height=300, background=vector(0, 0, 0))
bars = gvbars(delta=dv, color=color.green, label='Number of particles')
bars.plot(0, 0)
loopVerboseCount = 1

# Consts
dt = 5e-5
fps = 1000
k = 1.380649e-23
temperature = 300
globalStart = False

# Running
createWalls()
run(False)