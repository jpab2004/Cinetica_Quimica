#===============================================================================================#
#                                             Setup                                             #
#===============================================================================================#
# Libraries
# Iterable annotations
from collections.abc import Iterable

# Function to generate pseudorandom numbers equally distributed
from random import uniform, random

# Library for combinations of particles without repetition
from itertools import combinations

# Curve fitting
from lmfit.models import Model
from lmfit import Parameters

# Function to load element data for the simulation
from pickle import load

# Array support
from numpy import array
import numpy as np

# Function to verify if variable is not numeric
from math import isnan

# Functoin to pause script for time
from time import sleep

# Importing of all VPython tools (graphic library)
from vpython import *



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
    'histogramIterator': 1,
    'neighbourIterator': 1,
    'concentrationIterator': 1,
    'iterations': 1
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

fitFunction1X = lambda xs, a, b: [1/(a*x**b + 1) for x in xs]
fitFunctionInverted1X = lambda xs, a, b: [-1/(a*x**b + 1) + 1 for x in xs]

def fitFunctionE(x, amplitude, decay, independent):
    return amplitude*np.exp(decay*np.longdouble(x)) + independent



#===============================================================================================#
#                                        Velocities Graph                                       #
#===============================================================================================#
def generateTheoryCurveElement(e:int, nParticlesTheory:int) -> None:
    '''Generate the theory curve for a given number of particles of specific element.

    Args:
        e: int, element atomic number;
        nParticlesTheory: int, number of particles for the specified element.
    '''
    global velGraph, language

    theory = gcurve(color=elements[e]['color'], label=elements[e][language], graph=velGraph)
    mass = generateMassElement(e)

    deltaV = 10
    for v in range(0, maxVel+deltaV, deltaV):
        alpha = (dv**2/deltaV) * nParticlesTheory
        first = (mass / (2*pi*k*temperature))**1.5
        second = exp(-.5*mass*v**2 / (k*temperature))*v**2

        value = (alpha * first * second) / nParticles
        theory.plot(v, value)

    return



def generateTheoryCurveMolecule(mol:str, nParticlesTheory:int) -> None:
    '''Generate the theory curve for a given number of particles of specific molecule.

    Args:
        mol: str, molecule;
        nParticlesTheory: int, number of particles for the specified molecule.
    '''
    global velGraph, language

    theory = gcurve(color=molecules[mol]['color'], label=molecules[mol][language], graph=velGraph)    
    mass = generateMassMolecule(mol)

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



def drawHist(particles:Iterable[list, array]) -> None:
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

        if boxCorner:
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



def getRadiiElement(e:int) -> float:
    '''Get the atomic radius of given element.

    Args:
        e: int, element atomic number.
    
    Returns:
        float, radius of the particle with given element.
    '''
    global radiiFunc

    if empiricalRadii:
        radii = elements[e]['radii-empirical']/1e2
    else:
        radii = elements[e]['radii-calculated']/1e2

    if isnan(radii): return radiiFunc(150)    

    return radiiFunc(radii)



def getRadiiMolecule(mol:str) -> float:
    '''Get the atomic radius of given molecule.

    Args:
        mol: str, molecule.
    
    Returns:
        float, radius of the particle with given element.
    '''
    global radiiFunc

    if empiricalRadii:
        radii = molecules[mol]['radii-empirical']/1e2
    else:
        radii = molecules[mol]['radii-calculated']/1e2

    if isnan(radii): return radiiFunc(150)    

    return radiiFunc(radii)



def generateMassElement(e:int) -> float:
    '''Generates a mass for the particle of given element.

    Args:
        e: int, element atomic number.

    Returns:
        float, the mass of the particle of given element.
    '''
    return elements[e]['mass']*1E-3/6e23



def generateMassMolecule(mol:str) -> float:
    '''Generates a mass for the particle of given molecule.

    Args:
        mol: str, molecule.

    Returns:
        float, the mass of the particle of given molecule.
    '''
    return molecules[mol]['mass']*1E-3/6e23



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



def generateElement(e:int, d3:bool=True, overwriteTrail:bool=False, overwriteEmission:bool=False) -> sphere:
    '''Generate a particle (VPython Sphere object) for the simulation

    Args:
        e: int, element atomic number;
        d3: bool, True if the simulation is 3-Dimensional, false else.
    
    Returns:
        VPython Sphere object, the particle for the simulation.
    '''
    global makeEmission, makeTrails, retainTrail, elements

    particleRadius = getRadiiElement(e)
    particleMass = generateMassElement(e)
    particleColor = elements[e]['color']
    particleVelocity = generateVelocity(particleMass, d3)/particleMass
    particleTrail = ((overwriteTrail) or (makeTrails))
    particleEmission = ((overwriteEmission) or (makeEmission))

    still = True
    positions = list(map(getPositionAndRadius, globalParticles))
    while (still):
        still = False
        particlePosition = generatePosition(d3)
        for pos, rad in positions:
            if (mag(pos - particlePosition)) <= (rad + particleRadius):
                still = True
                break

    if prettySpheres:
        particle = sphere(pos=particlePosition, radius=particleRadius, color=particleColor, make_trail=particleTrail, retain=retainTrail)
    else:
        particle = simple_sphere(pos=particlePosition, radius=particleRadius, color=particleColor, makeTrails=particleTrail, retain=retainTrail)
    particle.v = particleVelocity
    particle.m = particleMass
    particle.emissive = particleEmission
    particle.type = e

    if neighbourImplementation:
        particle.neighbours = []
        particle.neighbourShell = particle.radius + mag(neighbourShellBuffer*dt*loopNeighboursCount*particle.v)
        particle.lastUpdatePosDistance = 0

    return particle



def generateMolecule(mol:str, d3:bool=True, overwriteTrail:bool=False, overwriteEmission:bool=False) -> sphere:
    '''Generate a particle (VPython Sphere object) for the simulation

    Args:
        mol: str, molecule;
        d3: bool, True if the simulation is 3-Dimensional, false else.
    
    Returns:
        VPython Sphere object, the particle for the simulation.
    '''
    global makeEmission, makeTrails, retainTrail, molecules

    particleRadius = getRadiiMolecule(mol)
    particleMass = generateMassMolecule(mol)
    particleColor = molecules[mol]['color']
    particleVelocity = generateVelocity(particleMass, d3)/particleMass
    particleTrail = ((overwriteTrail) or (makeTrails))
    particleEmission = ((overwriteEmission) or (makeEmission))

    still = True
    positions = list(map(getPositionAndRadius, globalParticles))
    while (still):
        still = False
        particlePosition = generatePosition(d3)
        for pos, rad in positions:
            if (mag(pos - particlePosition)) <= (rad + particleRadius):
                still = True
                break

    if prettySpheres:
        particle = sphere(pos=particlePosition, radius=particleRadius, color=particleColor, make_trail=particleTrail, retain=retainTrail)
    else:
        particle = simple_sphere(pos=particlePosition, radius=particleRadius, color=particleColor, make_trail=particleTrail, retain=retainTrail)
    particle.v = particleVelocity
    particle.m = particleMass
    particle.emissive = particleEmission
    particle.type = mol

    if neighbourImplementation:
        particle.neighbours = []
        particle.neighbourShell = particle.radius + mag(neighbourShellBuffer*dt*loopNeighboursCount*particle.v)
        particle.lastUpdatePosDistance = 0

    return particle



def generateFollowParticles(d3:bool=True):
    '''Generates the particles that can be follower via emission and trails.
    
    Args:
        d3: bool, True if the simulation is 3-Dimensional, false else.
    '''
    global globalParticles, followParticleTrail, followParticleEmission

    for element in elementsToSimulate:
        if elementsToSimulate[element] <= 0: continue
        elementsToSimulate[element] -= 1
        particle = generateElement(element, d3, overwriteTrail=followParticleTrail, overwriteEmission=followParticleEmission)
        globalParticles.append(particle)

    for mol in moleculesToSimulate:
        if moleculesToSimulate[mol] <= 0: continue
        moleculesToSimulate[mol] -= 1
        particle = generateMolecule(mol, d3, overwriteTrail=followParticleTrail, overwriteEmission=followParticleEmission)
        globalParticles.append(particle)

    return



def generateParticlesAndCurves(d3):
    '''Generate all particles and curves of the simulation.
    
    Args:
        d3: bool, True if the simulation is 3-Dimensional, false else.
    '''
    global elements, elementsToSimulate, molecules, moleculesToSimulate, globalTheoryCurve, globalParticles, showHistogram
    global followParticleTrail, followParticleEmission

    logicalFollow = ((followParticleTrail) or (followParticleEmission))
    if logicalFollow:
        generateFollowParticles(d3)

    for element, count in elementsToSimulate.items():
        if ((showHistogram) and (globalTheoryCurve)): generateTheoryCurveElement(element, (count + int(((logicalFollow) and (count > 0)))))
        for _ in range(count):
            particle = generateElement(element, d3)
            globalParticles.append(particle)

    for mol, count in moleculesToSimulate.items():
        if ((showHistogram) and (globalTheoryCurve)): generateTheoryCurveMolecule(mol, (count + int(((logicalFollow) and (count > 0)))))
        for _ in range(count):
            particle = generateMolecule(mol, d3)
            globalParticles.append(particle)
    
    return



def generateConcentrations():
    global concentrations, nParticles, conGraph, language
    global elements, elementsToSimulate, molecules, moleculesToSimulate

    for e, count in elementsToSimulate.items():
        curve = gcurve(graph=conGraph, color=elements[e]['color'], label=elements[e][language])
        concentrations[e] = [count, curve]
    for mol, count in moleculesToSimulate.items():
        curve = gcurve(graph=conGraph, color=molecules[mol]['color'], label=molecules[mol][language])
        concentrations[mol] = [count, curve]



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



def react(mol:str, p1:sphere, p2:sphere, toKill:Iterable[list, array]) -> Iterable[list, array]:
    '''Reacts two particles to make a molecule.
    
    Args:
        mol: String, molecule that the reaction is turning into;
        p1: VPython Sphere object, particle 1 of the reaction;
        p2: VPython Sphere object, particle 2 of the reaction;
        toKill: VPython Sphere object iterator, contains all particle objects to be killed.
    '''
    global molecules, concentrations

    concentrations[p1.type][0] -= 1
    concentrations[p2.type][0] -= 1
    concentrations[mol][0] += 1

    v1, v2 = p1.v, p2.v
    m1, m2 = p1.m, p2.m
    r1, r2 = p1.pos, p2.pos

    m_ = generateMassMolecule(mol)
    radius_ = getRadiiMolecule(mol)
    r_ = (r1 + r2) / 2
    v_ = (v1*m1 + v2*m2) / m_

    p1.type = mol
    p2.type = None
    p1.color = molecules[mol]['color']
    p1.trail_color = molecules[mol]['color']

    p1.m = m_
    p1.v = v_
    p1.pos = r_
    p1.radius = radius_

    toKill.append(p2)
    return toKill



def colision(p1, p2, toKill):
    '''Makes everything necessary for the collision.
    
    Args:
        p1: VPython Sphere object, particle 1 of the collision;
        p2: VPython Sphere object, particle 2 of the collision;
        toKill: VPython Sphere object iterator, contains all particle objects to be killed.
    '''

    if (mag(p1.pos - p2.pos) <= p1.radius + p2.radius) and (mag( (p1.pos + p1.v*dt) - (p2.pos + p2.v*dt) ) < mag(p1.pos - p2.pos)):
        for mol in moleculesToSimulate:
            if molecules[mol]['delay'] >= manager['iterations']:
                continue

            for type1, type2, chance in molecules[mol]['reagents-chance']:
                if ((chance >= random()) and (((type1 == p1.type) and (type2 == p2.type)) or ((type1 == p2.type) and (type2 == p1.type)))):
                    react(mol, p1, p2, toKill)
                    continue
    
        p1.v, p2.v = newVelocity(p1, p2)


def collisionDetection(particles:Iterable[list, array]) -> None:
    '''Detects and applies collisions for all particles.

    Args:
        particles: VPython Sphere object iterator, contains all particle objects in the simulation.
    '''
    global nParticles, moleculesToSimulate
    toKill = []

    if neighbourImplementation:
        for p1 in particles:
            for p2 in p1.neighbours:
                colision(p1, p2, toKill)
    else:
        for p1, p2 in combinations(particles, 2):
            if (mag(p1.pos - p2.pos) <= p1.radius + p2.radius) and (mag( (p1.pos + p1.v*dt) - (p2.pos + p2.v*dt) ) < mag(p1.pos - p2.pos)):
                colision(p1, p2, toKill)
    
    toKill = set(toKill)
    for p in toKill:
        p.clear_trail()
        particles.remove(p)
        nParticles -= 1
        p.visible = False
        del p
    
    return



def updateConcentrations():
    '''Updates the concentrations of all elements.'''
    global concentrations, manager

    for count, curve in concentrations.values():
        curve.plot(manager['iterations'], count/nParticles)

    return



def updateNeighbours(p1:sphere, particles:Iterable[list, array]) -> None:
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



def updateNeighboursAllParticles(particles:Iterable[list, array]) -> None:
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



def createFitCurve(params):
    '''Creates a fit graph of the simulation.
    
    Args:
        params: list, the parameters of the fit curve.
    '''
    global fittingGraph, language, fitCurveStopDelay, globalGdotRadius

    particleType, result, (xs, ys) = params

    if particleType in elements:
        fittedCurve = gcurve(color=elements[particleType]['color'], label=elements[particleType][language], graph=fittingGraph)
        dotsCurve = gdots(color=elements[particleType]['color'], graph=fittingGraph, radius=globalGdotRadius)
    elif particleType in molecules:
        fittedCurve = gcurve(color=molecules[particleType]['color'], label=molecules[particleType][language], graph=fittingGraph)
        dotsCurve = gdots(color=molecules[particleType]['color'], graph=fittingGraph, radius=globalGdotRadius)

    fittedCurve.plot([[x, y] for x, y in zip(xs, result.best_fit)])
    dotsCurve.plot([[x, y] for (x, y) in zip(xs, ys)])



def fitGraph():
    '''Fits the graph of the concentrations of the simulation.'''
    global concentrations, particlesToFit

    for p, (curveType, initial) in particlesToFit.items():
        x = [i[0] for i in concentrations[p][-1].data]
        y = [i[1] for i in concentrations[p][-1].data]

        model = Model(curveType)
        params = Parameters()

        for name, (value, vary, minv, maxv, expr, brute_step) in initial.items():
            params.add(name=name, value=value, vary=vary, min=minv, max=maxv, expr=expr, brute_step=brute_step)

        result = model.fit(y, params, x=x)

        print(result.fit_report())
        params = [p, result, [x, y]]
        createFitCurve(params)



#===============================================================================================#
#                                          3D Functions                                         #
#===============================================================================================#
def step3D(particles:Iterable[list, array]) -> None:
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
def step2D(particles:Iterable[list, array]) -> None:
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
def stepSimulation():
    '''Make 1 step in the simulation globaly.'''
    global showHistogram, showConcentrations, manager, globalFitCurveStop, fitCurveStopDelay, particlesToFit

    if ((neighbourImplementation) and (globalUpdateNeighbour) and (manager['neighbourIterator'] >= loopNeighboursCount)):
        updateNeighboursAllParticles(globalParticles)
        manager['neighbourIterator'] = 1

    rate(fps)
    manager['runFunction'](globalParticles)
    collisionDetection(globalParticles)

    if ((showConcentrations) and (manager['concentrationIterator'])) >= loopConcentrationVerboseCount:
        updateConcentrations()
        manager['concentrationIterator'] = 1
    if ((showHistogram) and (manager['histogramIterator'])) >= loopHistogramVerboseCount:
        drawHist(globalParticles)
        manager['histogramIterator'] = 1
    if ((globalFitCurveStop) and (manager['iterations'] >= fitCurveStopDelay)):
        fitGraph()
        globalFitCurveStop = False

    manager['histogramIterator'] += 1
    manager['concentrationIterator'] += 1
    if (neighbourImplementation): manager['neighbourIterator'] += 1

    manager['iterations'] += 1

    return



def run(d3:bool=True) -> None:
    '''Runs the global manager for the simulation.

    Observation:
        Creates all particle objects, prepares the initialization of the simulation
        and runs all necessary functions for the simulation step.

    Args:
        d3: bool, True if the simulation is 3-Dimensional, false else.
    '''
    global globalStart, startSimulationButton, manager

    if d3: manager['runFunction'] = step3D
    else: manager['runFunction'] = step2D

    generateParticlesAndCurves(d3)
    generateConcentrations()

    startSimulationButton.disabled = False

    while(not globalStart):
        rate(100)
        continue

    while(True):
        while(globalStart):
            while(not paused):
                stepSimulation()
        sleep(.1)



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



def stepSimulationNTimes():
    '''Make n step in the simulation globaly.'''
    global globalParticles, manager, stepSimulationButton, resumeSimulationButton

    stepSimulationButton.disabled = True
    resumeSimulationButton.disabled = True

    for _ in range(manager['numberOfSteps']):
        rate(50)
        stepSimulation()

    stepSimulationButton.disabled = False
    resumeSimulationButton.disabled = False

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
stepSimulationButton = button(pos=scene.caption_anchor, text='Rodar passos', bind=stepSimulationNTimes, disabled=True, left=50)

# Creates the slider to change the number of manual steps
numberOfStepsSlider = slider(min=1, max=500, value=250, step=1, length=220, bind=setNumberOfSteps, disabled=True, left=20)
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
temperature = 300
# Determine the language of the simulation
language = 'name-pt'
# Global start variable (global manager)
globalStart = False
# Global pause variable (global manager)
paused = False
# Neighbour otimization
neighbourImplementation = True
# Define if neighbours are update together or separetly (global manager)
# DO NOT CHANGE!!! NOT IMPLEMENTED CORRECTLY! WILL MAKE SIMULATION RUN SLOWER (MUCH SLOWER)!
globalUpdateNeighbour = True
# Defines if the histogram and velocities graph is created and shown
showHistogram = True
# Defines if the concentrations graph is created and shown
showConcentrations = generateTheoryCurveElement



# Wall variables
# Size of each wall (Angstrom)
side = 9
# Thickness of each wall
thickness = .5



# Emission and trail options
# Makes particles have trails and their retain iterations
makeTrails = False
retainTrail = 10
# Make particles emit light
makeEmission = False
# Make only 1 particle have a tral
followParticleTrail = True
# Make only 1 particle emit light
followParticleEmission = True



# Prettier
# Makes the simulation have prettier graphics (WIP)
boxCorner = True
# Make Vpython use prettier spheres
prettySpheres = True
# Defines if graphs use the fast implementation or not
globalFastGraph = False
# Define the graph background and foreground color
globalGraphBackgroundColor = HEX2VEC('#BBBBBB')
globalGraphForegroundColor = HEX2VEC('#000000')
# Citing the theory curve to be imprecise
scene.append_to_caption("<br><b>Obs.:</b> A curva teórica apresentada de Boltzmann esta relacionada ao estado inicial da simulação")
# Creating a line between graphs and controls
scene.append_to_caption("<br><hr><br>")



# Particle variables
# Particle list initialization
globalParticles = []
# Buffer for the radius size (graphics)
radiiFunc = lambda x: x/5 + .2
# Buffer for generating the initial position of particles
positionBuffer = .8*side
# Bool for defining if particles start randomly scattered or on the center
randomPosition = True
# Bool for defining use of empirical radii or calculated radii
empiricalRadii = True
# Amount of loops to update the list of neighbours of each particle (75 seens good for most pruposes) (DOES NOT WORK WITH LOW NUMBERS)
loopNeighboursCount = 75
# Size of the neighbour shell
neighbourShellBuffer = 1.5
# Size of the outer neighbour shell
neighbourMaxShellBuffer = 2.2



# Elements and Molecules
# List of element to simulate (atomic number)
elementsToSimulate = {1: 300}
# List of molecules to simulate
moleculesToSimulate = {'H2': 100}
# Total number of particles
nParticles = sum(elementsToSimulate.values()) + sum(moleculesToSimulate.values())



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
                ymax=1, ytitle='Densidade de Probabilidade', fast=globalFastGraph, width=velGraphWidth, height=velGraphHeight,
                align='left', background=globalGraphBackgroundColor, foreground=globalGraphForegroundColor)
# Velocities graph initialization
velBars = gvbars(delta=dv, color=HEX2VEC('#ff0000'), label='Número de Partículas', graph=velGraph)
# Amount of loops to update the graph (optimization)
loopHistogramVerboseCount = 5



# Reaction
# Lit of all concentrations
concentrations = {}
# Width and height of the concentration graph on the canvas
conGraphWidth, conGraphHeight = 850, 300
# Create a scrollable graph or not
scrollableConGraph = True
# Concentration graph creation
if scrollableConGraph:
    conGraph = graph(title='Concentração dos elementos/moléculas na simulação', xtitle='Tempo (iterações)',
    ytitle='Concentração', fast=globalFastGraph, width=conGraphWidth, height=conGraphHeight, align='left', ymin=0, ymax=1,
    background=globalGraphBackgroundColor, foreground=globalGraphForegroundColor, scroll=True, xmin=0, xmax=4*fps)
else:
    conGraph = graph(title='Concentração dos elementos/moléculas na simulação', xtitle='Tempo (iterações)',
    ytitle='Concentração', fast=globalFastGraph, width=conGraphWidth, height=conGraphHeight, align='left', ymin=0, ymax=1,
    background=globalGraphBackgroundColor, foreground=globalGraphForegroundColor)
# Amount of loops to update the graph (optimization)
loopConcentrationVerboseCount = 25



# Curve fitting
# Define the stop simulation iteration
globalFitCurveStop = True
fitCurveStopDelay = 2e4
# Defining the global size of plots
globalGdotRadius = 1
# Define fitting buffer
globalFitBuffer = 1e2
# Defining the particles to fit
particlesToFit = {
    1: (fitFunctionE, {
        'amplitude': [-1, True, -globalFitBuffer, globalFitBuffer, None, None],
        'decay': [-.004, True, -globalFitBuffer, 0, None, None],
        'independent': [1, True, -globalFitBuffer, globalFitBuffer, None, None]
    }),

    'H2': (fitFunctionE, {
        'amplitude':[1, True, -globalFitBuffer, globalFitBuffer, None, None],
        'decay': [-.004, True, -globalFitBuffer, 0, None, None],
        'independent': [0, True, -globalFitBuffer, globalFitBuffer, None, None]
    })
}
# Creating the graph
if globalFitCurveStop:
    fittingGraph = graph(title='Fit das concentrações na simulação', xtitle='Tempo (iterações)',
    ytitle='Concentração', fast=globalFastGraph, width=conGraphWidth, height=conGraphHeight, align='left', ymin=0, ymax=1,
    background=globalGraphBackgroundColor, foreground=globalGraphForegroundColor, xmin=0, xmax=fitCurveStopDelay)



# Running
# Creating the walls of the simulation
createWalls(False)
# Running the simulation
run(False)