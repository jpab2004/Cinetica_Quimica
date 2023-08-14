#===============================================================================================#
#                                             Setup                                             #
#===============================================================================================#
# Libraries
from vpython import scene, rate, box, sphere, cylinder, vector, button, canvas, dot, mag2
from pickle import load, HIGHEST_PROTOCOL
from math import cos, sin, radians, sqrt
from random import random, uniform
from itertools import combinations
from math import exp, isnan

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

def startSimulation():
    global globalStart
    globalStart = True
    startSimulationButton.disabled = True

    return

def euclidianDistance(p1, p2):
    return sqrt((p1.pos.x - p2.pos.x)**2 + (p1.pos.y - p2.pos.y)**2 + (p1.pos.z - p2.pos.z)**2)

def newVelocity(p1, p2):
    v1, v2 = p1.v, p2.v
    m1, m2 = p1.m, p2.m
    r1, r2 = p1.pos, p2.pos

    v1_ = v1 - ((2*m2) / (m1 + m2)) * (dot(v1 - v2, r1 - r2) / mag2(r1 - r2)) * (r1 - r2)
    v2_ = v2 - ((2*m1) / (m2 + m1)) * (dot(v2 - v1, r2 - r1) / mag2(r2 - r1)) * (r2 - r1)

    return v1_, v2_

def collision(iterator):
    for p1, p2 in combinations(iterator, 2):
        if (euclidianDistance(p1, p2) <= p1.radius + p2.radius):
            p1.v, p2.v = newVelocity(p1, p2)
    
    return

def getRadii(e):
    f = lambda x: exp(x/500) - .7

    if empiricalRadii:
        radii = elements[e]['radii-empirical']
    else:
        radii = elements[e]['radii-calculated']

    if isnan(radii): return f(150)    

    return f(radii)



#===============================================================================================#
#                                          3D Functions                                         #
#===============================================================================================#
def generate3DVelocity():
    if normalizedVelocity:
        psi = radians(uniform(0, 360))
        theta = radians(uniform(0, 360))

        x = normalizedVelocity * cos(theta) * sin(psi)
        y = normalizedVelocity * sin(theta) * sin(psi)
        z = normalizedVelocity * cos(psi)

        return vector(x, y, z)
    else:
        return vector.random()
    
def generate3DParticle(e):
    if randomPosition:
        position = vector(positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1))
    else:
        position = vector(0, 0, 0)

    particleRadius = getRadii(e)
    color = elements[e]['color']

    particle = sphere(pos=position, radius=particleRadius, color=color, make_trail=makeTrails, retain=10)
    
    particle.v = generate3DVelocity()
    particle.m = 2

    return particle

def step3D(iterator):
    for b in iterator:
        b.pos += b.v*dt
        wallCollision = side - .5*thickness - b.radius

        if not (wallCollision > b.pos.x > -wallCollision):
            b.v.x *= -1
        if not (wallCollision > b.pos.y > -wallCollision):
            b.v.y *= -1
        if not (wallCollision > b.pos.z > -wallCollision):
            b.v.z *= -1

    return



#===============================================================================================#
#                                          2D Functions                                         #
#===============================================================================================#
def generate2DVelocity():
    if normalizedVelocity:
        theta = radians(uniform(0, 360))

        x = normalizedVelocity * cos(theta)
        y = normalizedVelocity * sin(theta)

        return vector(x, y, 0)
    else:
        return vector(uniform(-1, 1), uniform(-1, 1), 0)

def generate2DParticle(e):
    if randomPosition:
        position = vector(positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1), 0)
    else:
        position = vector(0, 0, 0)

    particleRadius = getRadii(e)
    color = elements[e]['color']

    particle = sphere(pos=position, radius=particleRadius, color=color, make_trail=makeTrails, retain=10)
    
    particle.v = generate2DVelocity()
    particle.m = 2

    return particle

def step2D(iterator):
    for b in iterator:
        b.pos += (b.v/b.m)*dt
        wallCollision = side - .5*thickness - b.radius

        if not (wallCollision > b.pos.x > -wallCollision):
            b.v.x *= -1
            b.pos += (b.v/b.m)*dt
        if not (wallCollision > b.pos.y > -wallCollision):
            b.v.y *= -1
            b.pos += (b.v/b.m)*dt
    
    return



#===============================================================================================#
#                                       Running Functions                                       #
#===============================================================================================#
def run3D():
    particles = []
    for element, count in zip(elementsToSimulate, elementsCount):
        for _ in range(count):
            particle = generate3DParticle(element)
            particles.append(particle)

    while(not globalStart):
        rate(15)
        continue

    while(globalStart):
        rate(fps)
        step3D(particles)
        collision(particles)

    return

def run2D():
    particles = []
    for element, count in zip(elementsToSimulate, elementsCount):
        for _ in range(count):
            particle = generate2DParticle(element)
            particles.append(particle)

    while(not globalStart):
        rate(15)
        continue

    while(globalStart):
        rate(fps)
        step2D(particles)
        collision(particles)
    
    return



#===============================================================================================#
#                                          Controllers                                          #
#===============================================================================================#
startSimulationButton = button(pos=scene.caption_anchor, text='Start simulation', bind=startSimulation)



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
radiiBuff = .01
normalizedVelocity = 5
positionBuffer = 8
randomPosition = True
empiricalRadii = True

# Elements
elementsToSimulate = [1, 8]
elementsCount = [20, 20]

# Consts
fps = 20*sum(elementsCount)
dt = .01
globalStart = False

# Running
createWalls()
run3D()