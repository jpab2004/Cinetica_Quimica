#===============================================================================================#
#                                             Setup                                             #
#===============================================================================================#
# Libraries
# Iterable annotations
from collections.abc import Iterable

# Library for combinations of particles without repetition
from itertools import combinations

# Function to generate pseudorandom numbers equally distributed
from random import uniform, random

# Function to load element data for the simulation
from pickle import load

# Function to verify if variable is not numeric
from math import isnan

# Functoin to pause script for time
from time import sleep

# Importing of all VPython tools (graphic library)
from vpython import *

# Array support
import numpy



# Preparing Scene
# Deleting the pre-initialized scene
scene.delete()
# Defining the buffer for the dimensions
sceneBuffer = 1
# Defining the scene width and height
sceneWidth = 900 * sceneBuffer
sceneHeight = 900 * sceneBuffer
# Creating the canvas of the simulation
scene = canvas(width=sceneHeight, height=sceneHeight, align='left')
# Setting the background color
scene.background = vector(0, 0, 0)
# Appending the processing information to the caption of the simulation
scene.append_to_caption("<div id='fps'/>")



# GLOBAL MANAGER
# Global simulation manager! DON'T CHANGE
manager = {
    'i': 1,
    'j': 1
}



# Elements
# Importing the elements data
file1 = 'Assets/elementsHashTable.pickle'
with open(file1, 'rb') as f:
    elements = load(f)

file2 = 'Assets/moleculesHashTable.pickle'
with open(file2, 'rb') as f:
    molecules = load(f)



#===============================================================================================#
#                                         Aux Functions                                         #
#===============================================================================================#
def RGB2VEC(r, g, b): 
    return vector(r/255, g/255, b/255)

def HEX2VEC(hex):
    hex = hex.lstrip('#')
    r, g, b = tuple(int(hex[i:i+2], 16) for i in (0, 2, 4))
    return RGB2VEC(r, g, b)



#===============================================================================================#
#                                        Velocities Graph                                       #
#===============================================================================================#
def generateTheoryCurve(e:int, nParticlesTheory:int) -> None:
    '''Generate the theory curve for a given number of particles of specific element.

    Args:
        e: int, element atomic number;
        nParticlesTheory: int, number of particles for the specified element.
    '''
    global velGraph
    
    theory = gcurve(color=elements[e]['color'], label=elements[e]['name-pt'], graph=velGraph)
    mass = generateMass(e)

    deltaV = 10
    for v in range(0, maxVel+deltaV, deltaV):
        alpha = (dv**2/deltaV) * nParticlesTheory
        first = (mass / (2*pi*k*temperature))**1.5
        second = exp(-.5*mass*v**2 / (k*temperature))*v**2

        value = (alpha * first * second) / nParticles
        theory.plot(v, value)

    return



def getVelocity(p:sphere) -> float:
    '''Calculates magnitute of velocity of the particle.

    Args:
        p: VPpython Sphere object, the particle.
    
    Returns:
        float, magnitute of the velocity.
    '''
    return mag(p.v)



def getHist(v:float) -> int:
    '''Calculates the bin for the velocity based on deltaV.

    Args:
        v: float, velocity to calculate the bin.

    Returns:
        int, bin index of the velocity.
    '''
    return int(v/dv)



def drawHist(particles:Iterable[list, numpy.array]) -> None:
    '''Generates and plots the histogram of velocities of the simulation.

    Args:
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    global velBars

    histData = {}
    vels = list(map(getVelocity, particles))
    for i in list(map(getHist, vels)):
        try:
            histData[i] += 1
        except:
            histData[i] = 0

    histData = [[v*dv + .5*dv, histData[v]/nParticles] if v in histData else [v*dv + .5*dv, 0] for v in range(max(histData))]
    velBars.data = histData

    return



#===============================================================================================#
#                                        Common Functions                                       #
#===============================================================================================#
def RGB2VEC(r:float, g:float, b:float) -> vector:
    '''Generates a VPython Vector object from RGB information.

    Args:
        r: float, red channel for the color;
        g: float, green channel for the color;
        b: float, blue channel for the color.

    Returns:
        VPython Vector, vector representing the color.
    '''
    return vector(r/255, g/255, b/255)



def createWalls(solidWalls:bool) -> None:
    '''Generate the walls of the simulation.

    Args:
        solidWalls: bool, defines if the walls will be solid or wireframes.
    '''
    if solidWalls:
        s2 = 2 * side - thickness
        s3 = 2 * side + thickness

        wallE = box(pos=vector(side, 0, 0), size=vector(thickness, s2, s2), color=vector(1, 0, 0))
        wallW = box(pos=vector(-side, 0, 0), size=vector(thickness, s2, s2), color=vector(1, 0, 0))

        wallN = box(pos=vector(0, side, 0), size=vector(s3, thickness, s2), color=vector(0, 1, 0))
        wallS = box(pos=vector(0, -side, 0), size=vector(s3, thickness, s2), color=vector(0, 1, 0))
        
        wallC = box(pos=vector(0, 0, -side), size=vector(s3, s3, thickness), color=vector(0, 0, 1))
    else:
        edgeColor = HEX2VEC('#dbf54c')
        edgeRadius = side*.03
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



def getRadii(e:int) -> float:
    '''Get the atomic radius of given element.

    Args:
        e: int, element atomic number.
    
    Returns:
        float, radius of the particle with given element.
    '''
    f = lambda x: (exp(x/500) - .7) * radiiBuff

    if empiricalRadii:
        radii = elements[e]['radii-empirical']/1e2
    else:
        radii = elements[e]['radii-calculated']/1e2

    if isnan(radii): return f(150)    

    return f(radii)



def generateMass(e:int) -> float:
    '''Generates a mass for the particle with given element.

    Args:
        e: int, element atomic number.

    Returns:
        float, the mass of the particle with given element.
    '''
    return elements[e]['mass']*1E-3/6e23



def generateVelocity(particleMass:float, d3:bool=True) -> vector:
    '''Generates a velocity vector for particle initialization.

    Args:
        particleMass: float, mass of the particle to be initialized;
        d3: bool, True if the simulation is 3-Dimensional, false else.
    
    Returns:
        VPython vector object, vector with the velocity for the particle.
    '''
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



def getPositionAndRadius(p:sphere) -> tuple[vector, float]:
    '''Get the position and radius of particle.

    Args:
        p: VPython sphere object, particle to get position and radius.

    Returns:
        tuple, position and radius of particle.
    '''
    return (p.pos, p.radius)



def generatePosition(d3:bool=True) -> vector:
    '''Generates a postion vector for particle initialization.

    Args:
        d3: bool, True if the simulation is 3-Dimensional, false else.

    Returns:
        VPython Vector object, vector with position on the scene.
    '''
    if randomPosition:
        if d3: return positionBuffer*vector.random()
        else: return vector(positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1), 0)
    else: return vector(0, 0, 0)



def generateParticle(e:int, d3:bool=True) -> sphere:
    '''Generate a particle (VPython Sphere object) for the simulation

    Args:
        e: int, element atomic number;
        d3: bool, True if the simulation is 3-Dimensional, false else.
    
    Returns:
        VPython Sphere object, the particle for the simulation.
    '''
    global particleEmission

    particleRadius = getRadii(e)
    particleMass = generateMass(e)
    particleColor = elements[e]['color']
    particleVelocity = generateVelocity(particleMass, d3)/particleMass

    still = True
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
    particle.emissive = particleEmission
    particle.element = e
    particle.mol = None

    if neighbourImplementation:
        particle.neighbours = []
        particle.neighbourShell = particle.radius + mag(neighbourShellBuffer*dt*loopNeighboursCount*particle.v)
        particle.lastUpdatePosDistance = 0

    return particle



def newVelocity(p1:sphere, p2:sphere) -> tuple[vector, vector]:
    '''Calculates and returns the new velocity for 2 particles colliding.

    Args:
        p1: VPython Sphere object, particle 1 of the collision;
        p2: VPython Sphere object, particle 2 of the collision.

    Returns:
        tuple, particles velocities for p1 and p2 sequentially.
    '''
    v1, v2 = p1.v, p2.v
    m1, m2 = p1.m, p2.m
    r1, r2 = p1.pos, p2.pos

    v1_ = v1 - ((2*m2) / (m1 + m2)) * (dot(v1 - v2, r1 - r2) / mag2(r1 - r2)) * (r1 - r2)
    v2_ = v2 - ((2*m1) / (m2 + m1)) * (dot(v2 - v1, r2 - r1) / mag2(r2 - r1)) * (r2 - r1)

    return (v1_, v2_)



def react(mol:str, p1:sphere, p2:sphere, toKill:Iterable[list, numpy.array]) -> Iterable[list, numpy.array]:
    '''Reacts two particles to make a molecule.
    
    Args:
        mol: String, molecule that the reaction is turning into;
        p1: VPython Sphere object, particle 1 of the reaction;
        p2: VPython Sphere object, particle 2 of the reaction;
        toKill: VPython Sphere object iterator, contains all particle objects to be killed.
    '''
    v1, v2 = p1.v, p2.v
    m1, m2 = p1.m, p2.m
    r1, r2 = p1.pos, p2.pos
    radius1, radius2 = p1.radius, p2.radius

    m_ = m1 + m2
    radius_ = (radius1 + radius2) / 1.5
    r_ = (r1 + r2) / 2
    v_ = (v1*m1 + v2*m2) / m_

    p1.element = None
    p1.mol = mol
    p1.color = molecules[mol]['color']
    p1.trail_color = molecules[mol]['color']

    p1.m = m_
    p1.v = v_
    p1.pos = r_
    p1.radius = radius_

    toKill.append(p2)
    return toKill



def collision(particles:Iterable[list, numpy.array]) -> None:
    '''Detects and applies collisions for all particles.

    Args:
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    global nParticles
    toKill = []

    if neighbourImplementation:
        for p1 in particles:
            for p2 in p1.neighbours:
                if (mag(p1.pos - p2.pos) <= p1.radius + p2.radius) and (mag( (p1.pos + p1.v*dt) - (p2.pos + p2.v*dt) ) < mag(p1.pos - p2.pos)):
                    if ((H2ReactionChance >= random()) and (p1.element == p2.element == 1)):
                        toKill = react('H2', p1, p2, toKill)
                    else:
                        p1.v, p2.v = newVelocity(p1, p2)
    else:
        for p1, p2 in combinations(particles, 2):
            if (mag(p1.pos - p2.pos) <= p1.radius + p2.radius) and (mag( (p1.pos + p1.v*dt) - (p2.pos + p2.v*dt) ) < mag(p1.pos - p2.pos)):
                if ((H2ReactionChance >= random()) and (p1.element == p2.element == 1)):
                    toKill = react('H2', p1, p2, toKill)
                else:
                    p1.v, p2.v = newVelocity(p1, p2)

    for p in toKill:
        p.visible = False
        p.pos = vector(100,100,0)
        p.vel = vector(0,0,0)
        p.clear_trail()
        particles.remove(p)
        nParticles -= 1
        del p
    
    return



def updateNeighbours(p1:sphere, particles:Iterable[list, numpy.array]) -> None:
    '''Updates the neighbours of 1 particle.
    
    Args:
        p1: VPython Sphere object, particle to update
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    setattr(p1, 'neighbours', [])
    p1.neighbourShell = p1.radius + mag(neighbourShellBuffer*dt*loopNeighboursCount*p1.v)
    p1.lastUpdatePosDistance = 0
    for p2 in particles:
        if p1 == p2:continue
        if (mag(p1.pos - p2.pos) <= (p1.radius + neighbourMaxShellBuffer*dt*loopNeighboursCount*max(mag(p1.v), mag(p2.v)))):
            p1.neighbours.append(p2)
    
    return



def updateNeighboursAllParticles(particles:Iterable[list, numpy.array]) -> None:
    '''Updates the neighbours of each particles.
    
    Args:
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    [setattr(p, 'neighbours', []) for p in particles]
    for p1, p2 in combinations(particles, 2):   
        if (mag(p1.pos - p2.pos) <= (p1.radius + neighbourMaxShellBuffer*dt*loopNeighboursCount*max(mag(p1.v), mag(p2.v)))):
            p1.neighbours.append(p2)
            p2.neighbours.append(p1)

    return



#===============================================================================================#
#                                          3D Functions                                         #
#===============================================================================================#
def step3D(particles:Iterable[list, numpy.array]) -> None:
    '''Function to calculate next step of the simulation for 3D simulations.

    Args:
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    for p in particles:
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

        if ((neighbourImplementation) and (not globalUpdateNeighbour) and (p.lastUpdatePosDistance >= p.neighbourShell)): updateNeighbours(p, particles)

    return



#===============================================================================================#
#                                          2D Functions                                         #
#===============================================================================================#
def step2D(particles:Iterable[list, numpy.array]) -> None:
    '''Function to calculate next step of the simulation for 2D simulations.

    Args:
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    for p in particles:
        if neighbourImplementation: lastPos = mag(p.pos)

        p.pos += p.v*dt
        wallCollision = side - .5*thickness - p.radius

        if neighbourImplementation: p.lastUpdatePosDistance += abs(lastPos - mag(p.pos))

        if not (wallCollision > p.pos.x > -wallCollision):
            p.v.x *= -1
            p.pos += p.v*dt
        if not (wallCollision > p.pos.y > -wallCollision):
            p.v.y *= -1
            p.pos += p.v*dt

        if ((neighbourImplementation) and (not globalUpdateNeighbour) and (p.lastUpdatePosDistance >= p.neighbourShell)): updateNeighbours(p, particles)
    
    return



#===============================================================================================#
#                                       Running Function                                        #
#===============================================================================================#
def run(d3:bool=True) -> None:
    '''Runs the global manager for the simulation.

    Observation:
        Creates all particle objects, prepares the initialization of the simulation
        and runs all necessary functions for the simulation step.

    Args:
        d3: bool, True if the simulation is 3-Dimensional, false else.
    '''
    global globalStart, startSimulationButton, particles, manager

    if d3: manager['runFunction'] = step3D
    else: manager['runFunction'] = step2D

    for element, count in zip(elementsToSimulate, elementsCount):
        if globalTheoryCurve: generateTheoryCurve(element, count)
        for _ in range(count):
            particle = generateParticle(element, d3)
            particles.append(particle)

    startSimulationButton.disabled = False

    while(not globalStart):
        rate(15)
        continue

    while(True):
        while(globalStart):
            while(not paused):
                if ((neighbourImplementation) and (globalUpdateNeighbour) and (manager['j'] >= loopNeighboursCount)): updateNeighboursAllParticles(particles); manager['j'] = 1

                rate(fps)
                manager['runFunction'](particles)
                collision(particles)

                if manager['i'] >= loopVerboseCount: drawHist(particles); manager['i'] = 1

                manager['i'] = manager['i'] + 1
                if (neighbourImplementation): manager['j'] = manager['j'] + 1
        sleep(.1)
    
    return



#===============================================================================================#
#                                          Controllers                                          #
#===============================================================================================#
def startSimulation():
    '''Starts the simulation globaly.'''
    global globalStart, startSimulationButton, pauseSimulationButton, numberOfStepsSlider

    globalStart = True
    startSimulationButton.disabled = True
    pauseSimulationButton.disabled = False
    numberOfStepsSlider.disabled = True

    return



def pauseSimulation():
    '''Pauses the simulation globaly.'''
    global paused, pauseSimulationButton, resumeSimulationButton
    global stepSimulationButton, numberOfStepsSlider
    
    paused = True
    pauseSimulationButton.disabled = True
    resumeSimulationButton.disabled = False
    stepSimulationButton.disabled = False
    numberOfStepsSlider.disabled = False

    return



def resumeSimulation():
    '''Resumes the simulation globaly.'''
    global paused, pauseSimulationButton, resumeSimulationButton
    global stepSimulationButton, numberOfStepsSlider

    paused = False
    pauseSimulationButton.disabled = False
    resumeSimulationButton.disabled = True
    stepSimulationButton.disabled = True
    numberOfStepsSlider.disabled = True

    return



def stepSimulation():
    '''Make 1 step in the simulation globaly.'''
    global particles, manager

    for _ in range(manager['numberOfSteps']):
        if ((neighbourImplementation) and (globalUpdateNeighbour) and (manager['j'] >= loopNeighboursCount)): updateNeighboursAllParticles(particles); manager['j'] = 1

        rate(fps)
        manager['runFunction'](particles)

        if manager['i'] >= loopVerboseCount: drawHist(particles); manager['i'] = 1
        manager['i'] = manager['i'] + 1
        if (neighbourImplementation): manager['j'] = manager['j'] + 1

    return



def setNumberOfSteps(s:slider) -> None:
    '''Set the number of steps of the manual simulation.'''
    global manager, numberOfStepsText

    x = s.value
    numberOfStepsText.text = f'Número de passos: {x}'
    manager['numberOfSteps'] = x

    return



# Creates the button to start the simulation
startSimulationButton = button(pos=scene.caption_anchor, text='Começar simulação', bind=startSimulation, disabled=True, left=50)
# Creates the button to pause the simulation
pauseSimulationButton = button(pos=scene.caption_anchor, text='Pausar simulação', bind=pauseSimulation, disabled=True, left=50)

# Appending a break of line to the caption
scene.append_to_caption("<br>")

# Creates the button to resume the simulation
resumeSimulationButton = button(pos=scene.caption_anchor, text='Resumir simulação', bind=resumeSimulation, disabled=True, left=50)
# Creates the button to make 1 step in the simulation
stepSimulationButton = button(pos=scene.caption_anchor, text='Rodar passos', bind=stepSimulation, disabled=True, left=50)

# Creates the slider to change the number of manual steps
numberOfStepsSlider = slider(min=1, max=100, value=20, step=1, length=220, bind=setNumberOfSteps, disabled=True, left=20)
numberOfStepsText = wtext(text=f'Número de passos: {numberOfStepsSlider.value}')
manager['numberOfSteps'] = numberOfStepsSlider.value



#===============================================================================================#
#                                           Simulation                                          #
#===============================================================================================#
# Consts
# Delta Time for steps on the simulation
dt = 5e-6
# FPS of the simulation (max available fps, can be less)
fps = 5000
# Boltzmann Constant for calculations
k = 1.380649e-23
# Temperature of the simulation
temperature = 100
# Global start variable (global manager)
globalStart = False
# Global pause variable (global manager)
paused = False
# Neighbour otimization
neighbourImplementation = True
# Define if neighbours are update together or separetly (global manager)
# DO NOT CHANGE!!! NOT IMPLEMENTED CORRECTLY! WILL MAKE SIMULATION RUN SLOWER (MUCH SLOWER)!
globalUpdateNeighbour = True



# Wall variables
# Size of each wall (Angstrom)
side = 20
# Thickness of each wall
thickness = .5



# Prettier
# Makes the simulation have prettier graphics (WIP)
overallPretty = True
# Makes particles have trails
makeTrails = False
# Make particles emit light
particleEmission = False
# Creating a line between graphs and controls
scene.append_to_caption("<br><hr><br>")



# Particle variables
# Particle list initialization
particles = []
# Buffer for the radius size (graphics)
radiiBuff = 1
# Buffer for generating the initial position of particles
positionBuffer = .9*side
# Bool for defining if particles start randomly scattered or on the center
randomPosition = True
# Bool for defining use of empirical radii or calculated radii
empiricalRadii = True
# Amount of loops to update the list of neighbours of each particle (DOES NOT WORK WITH LOW NUMBERS)
loopNeighboursCount = 75
# Size of the neighbour shell
neighbourShellBuffer = 1.5
# Size of the outer neighbour shell
neighbourMaxShellBuffer = 2.2



# Elements
# List of element to simulate (atomic number)
elementsToSimulate = [1]
# The amount of each element to simulate
elementsCount = [300]
# Total number of particles
nParticles = sum(elementsCount)



# Velocity graph variables
# Determine whether to generate the theory curve of each element or not
globalTheoryCurve = True
# Delta V for histogram binning
dv = 100
# Max velocity displayed on the graph
maxVel = 6000
# Width and height of the velocity graph on the canvas
velGraphWidth, velGraphHeight = 850, 300
# Initial histogram data for plotting
histData = [(v*dv + .5*dv, 0) for v in range(int(maxVel/dv))]
# Velocities graph creation
velGraph = graph(title='Velocidade das partículas na simulação', xtitle='Velocidade (m/s)', xmax=maxVel,
                ymax=1, ytitle='Densidade de Probabilidade', fast=False, width=velGraphWidth, height=velGraphHeight,
                align='left', background=vector(0, 0, 0), foreground=vector(0, 0, 0))
# Velocities graph initialization
velBars = gvbars(delta=dv, color=HEX2VEC('#00ff00'), label='Número de Partículas', graph=velGraph)
velBars.plot(0, 0)
# Amount of loops to update the graph (optimization)
loopVerboseCount = 5



# Reaction
# Chance for a reaction to happen
H2ReactionChance = .1
# Width and height of the concentration graph on the canvas
conGraphWidth, conGraphHeight = 850, 300
# Create a scrollable graph or not
scrollableConGraph = True
# Concentration graph creation
if scrollableConGraph:
    conGraph = graph(title='Concentração dos elementos/moléculas na simulação', xtitle='Tempo (iterações)',
    ytitle='Concentração', fast=False, width=conGraphWidth, height=conGraphHeight, align='left',
    background=vector(0, 0, 0), foreground=vector(0, 0, 0), scroll=True, xmin=0, xmax=10*fps)
else:
    conGraph = graph(title='Concentração dos elementos/moléculas na simulação', xtitle='Tempo (iterações)',
    ytitle='Concentração', fast=False, width=conGraphWidth, height=conGraphHeight, align='left',
    background=vector(0, 0, 0), foreground=vector(0, 0, 0))
g = gcurve(graph=conGraph, color=HEX2VEC('#ff0000'))
g.plot(0, 1)
g.plot(5*fps, 2)



# Running
# Creating the walls of the simulation
createWalls(False)
# Running the simulation
run(False)