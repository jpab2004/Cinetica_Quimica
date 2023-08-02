#===============================================================================================#
#                                             Setup                                             #
#===============================================================================================#
# Libraries
from vpython import scene, box, sphere, cylinder, arrow, vector, rate, dot
from random import random, uniform
from math import cos, sin, radians
from time import sleep



#===============================================================================================#
#                                        Common Functions                                       #
#===============================================================================================#
RGB2VEC = lambda r, g, b: vector(r/255, g/255, b/255)

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
    
def generate3DBall():
    if randomPosition:
        positionBuffer = 3
        position = vector(positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1))
    else:
        position = vector(0, 0, 0)

    if randomColor:
        colorBuffer = 1.5
        color = vector(colorBuffer*random(), colorBuffer*random(), colorBuffer*random())
    else:
        color = vector(.8, .8, .8)

    ball = sphere(pos=position, radius=ballRadius, color=color, make_trail=makeTrails, retain=10)
    
    ball.v = generate3DVelocity()
    ball.m = 2

    return ball

def step3D(iterator):
    for b in iterator:
        b.pos += (b.v/b.m)*dt

        if not (wallCollision > b.pos.x > -wallCollision):
            b.v.x *= -1
        if not (wallCollision > b.pos.y > -wallCollision):
            b.v.y *= -1
        if not (wallCollision > b.pos.z > -wallCollision):
            b.v.z *= -1

    return

def collision3D(iterator):
    for i, b1 in enumerate(iterator):
        for b2 in iterator[i+1:]:
            xCol = abs(b1.pos.x - b2.pos.x) < b1.radius + b2.radius
            yCol = abs(b1.pos.y - b2.pos.y) < b1.radius + b2.radius
            zCol = abs(b1.pos.z - b2.pos.z) < b1.radius + b2.radius
            if (xCol and yCol and zCol):
                print('Collision!')
    
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

def generate2DBall():
    if randomPosition:
        positionBuffer = 3
        position = vector(positionBuffer*uniform(-1, 1), positionBuffer*uniform(-1, 1), 0)
    else:
        position = vector(0, 0, 0)

    if randomColor:
        colorBuffer = 1.5
        color = vector(colorBuffer*random(), colorBuffer*random(), colorBuffer*random())
    else:
        color = vector(.8, .8, .8)

    ball = sphere(pos=position, radius=ballRadius, color=color, make_trail=makeTrails, retain=50)
    
    ball.v = generate2DVelocity()
    ball.m = 2

    ball.unitNormal = arrow(round=True, pos=ball.pos, axis=vector(0, 0, 0), length=0, color=RGB2VEC(252, 152, 3))
    ball.unitTangent = arrow(round=True, pos=ball.pos, axis=vector(0, 0, 0), length=0, color=RGB2VEC(244, 252, 3))
    ball.unitVelocity = arrow(round=True, pos=ball.pos, axis=ball.v, length=0, color=RGB2VEC(3, 252, 240))
    ball.unitNewVelocity = arrow(round=True, pos=ball.pos, axis=ball.v, length=0, color=RGB2VEC(3, 252, 240))

    return ball

def step2D(iterator):
    for b in iterator:
        b.pos += (b.v/b.m)*dt

        if not (wallCollision > b.pos.x > -wallCollision):
            b.v.x *= -1
        if not (wallCollision > b.pos.y > -wallCollision):
            b.v.y *= -1
    
    return

def collision2D(iterator):
    for i, b1 in enumerate(iterator):
        for b2 in iterator[i+1:]:
            xCol = abs(b1.pos.x - b2.pos.x) < b1.radius + b2.radius
            yCol = abs(b1.pos.y - b2.pos.y) < b1.radius + b2.radius
            if (xCol and yCol):
                print('Collision!')
                # b1.v = b1.unitNewVelocity.axis
                # b2.v = b2.unitNewVelocity.axis

def draw2D(iterator):
    lengthBuffer = 2
    for i, b1 in enumerate(iterator):
        for b2 in iterator[i+1:]:
            x = b2.pos.x - b1.pos.x
            y = b2.pos.y - b1.pos.y

            b1.unitNormal.pos = b1.pos
            b1.unitNormal.axis = vector(x, y, 0)
            b1.unitNormal.length = lengthBuffer*b1.radius

            b2.unitNormal.pos = b2.pos
            b2.unitNormal.axis = vector(-x, -y, 0)
            b2.unitNormal.length = lengthBuffer*b2.radius

            b1.unitTangent.pos = (b1.pos + b2.pos)/2
            b1.unitTangent.axis = vector(-y, x, 0)
            b1.unitTangent.length = lengthBuffer*b1.radius

            b2.unitTangent.pos = (b2.pos + b1.pos)/2
            b2.unitTangent.axis = vector(y, -x, 0)
            b2.unitTangent.length = lengthBuffer*b2.radius

            b1.unitVelocity.pos = b1.pos
            b1.unitVelocity.axis = b1.v
            b1.unitVelocity.length = lengthBuffer*b1.radius

            b2.unitVelocity.pos = b2.pos
            b2.unitVelocity.axis = b2.v
            b2.unitVelocity.length = lengthBuffer*b2.radius

            # v1n = dot(b1.v, b1.unitNormal.axis)
            # v1t = dot(b1.v, b1.unitTangent.axis)
            # v2n = dot(b2.v, b2.unitNormal.axis)
            # v2t = dot(b2.v, b2.unitTangent.axis)

            # v1n_ = (v1n * (b1.m - b2.m) + 2 * b2.m * v2n) / (b1.m + b2.m)
            # v1t_ = v1t
            # v2n_ = (v2n * (b2.m - b1.m) + 2 * b1.m * v1n) / (b2.m + b1.m)
            # v2t_ = v2t

            # v1n_ *= b1.unitNormal.axis
            # v1t_ *= b1.unitTangent.axis
            # v2n_ *= b2.unitNormal.axis
            # v2t_ *= b2.unitTangent.axis

            # v1_ = v1n_ + v1t_
            # v2_ = v2n_ + v2t_

            # b1.unitNewVelocity.pos = (b1.pos + b2.pos)/2
            # b1.unitNewVelocity.axis = v1_
            # b1.unitNewVelocity.length = lengthBuffer*b1.radius

            # b2.unitNewVelocity.pos = (b2.pos + b1.pos)/2
            # b2.unitNewVelocity.axis = v2_
            # b2.unitNewVelocity.length = lengthBuffer*b2.radius



#===============================================================================================#
#                                       Running Functions                                       #
#===============================================================================================#
def run3D():
    balls = []
    for _ in range(ballCount):
        ball = generate3DBall()
        balls.append(ball)

    play = True

    while(play):
        rate(fps)
        step3D(balls)
        collision3D(balls)

    return

def run2D():
    balls = []
    for _ in range(ballCount):
        ball = generate2DBall()
        balls.append(ball)

    play = True
    
    while(play):
        rate(fps)
        step2D(balls)
        draw2D(balls)
        collision2D(balls)
    
    return



#===============================================================================================#
#                                           Simulation                                          #
#===============================================================================================#
# Preparing Scene
sceneBuffer = .8
scene.width = 1920 * sceneBuffer
scene.height = 1080 * sceneBuffer
scene.background = vector(0, 0, 0)

# Wall variables
side = 10
thickness = .5

# Ball variables
ballCount = 2
ballRadius = 1
normalizedVelocity = .5

# Ball setup
randomPosition = True
wallCollision = side - .5*thickness - ballRadius

# Prettier
overallPretty = False
randomColor = False
makeTrails = False
solidWalls = False

# Consts
fps = 300
dt = .1

# Running
createWalls()
sleep(1)
run2D()