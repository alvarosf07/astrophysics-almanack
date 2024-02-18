# Step #1 - Creation of 3 satellites (LEO,MEO,GEO) orbiting the Earth

from __future__ import division, print_function
from vpython import *

# Right button drag or Ctrl-drag to rotate "camera" to view scene. This program only works with vpython
# To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
#  On a two-button mouse, middle is left + right.

scene=canvas(title="VPython Example")

G = 6.67e-11

# First, we define the characteristics of the earth (initial position vector
# radius, color, with/without trail, kind of trail, interval of refreshing)
earth = sphere(pos=vector(0,0,0), radius=6.4e6, color=color.white, make_trail=True, trail_type='points', interval=10, retain=20,texture={'file':textures.earth})
earth.mass = 5.98e24
earth.p = vector(0, 0, 0) * earth.mass

# We repeat the same process for the MEO satellite
satellite = sphere(pos=vector(0,(6.4e6+20e6),0), radius=1e6, color=color.yellow,make_trail=True, retain=200)
satellite.speed = vector((sqrt(G*earth.mass/(6.4e6+20e6))),0,0)

# We repeat the same process for the LEO satellite
satellite2 = sphere(pos=vector(0,(6.4e6+163e3),0), radius=1e6, color=color.red,make_trail=True, retain=20)
satellite2.speed = vector((sqrt(G*earth.mass/(6.4e6+163e3))),0,0)

# We repeat the same process for the GEO satellite
satellite3 = sphere(pos=vector(0,(6.4e6+36e6),0), radius=1e6, color=color.green, make_trail=True, retain=200)
satellite3.speed = vector((sqrt(G*earth.mass/(6.4e6+36e6))),0,0)

# We start the infinite loop that will compute the star movement.
# We first define the time interval (dt). Then the loop starts.
dt = 10
while True:
    rate(100)
    # Definition of r as the vector that indicates the difference of distances 
    r = earth.pos - satellite.pos
    r2 = earth.pos - satellite2.pos
    r3 = earth.pos - satellite3.pos
    
    # Calculation of the Force that the earth exerts over the satellite (Newton Equation)
    F = (G * earth.mass * norm(r)) / (mag2(r)) # Units: Force/mass
    F2 = (G * earth.mass * norm(r2)) / (mag2(r2)) # Units: Force/mass
    F3 = (G * earth.mass * norm(r3)) / (mag2(r3)) # Units: Force/mass
    
    # Calculation of the speed vector of the satellites at each moment (dt)
    satellite.speed = satellite.speed + F*dt 
    satellite2.speed = satellite2.speed + F2*dt
    satellite3.speed = satellite3.speed + F3*dt
    
    # Calculation of the position of the satellite at each moment (dt)
    satellite.pos = satellite.pos + (satellite.speed * dt)
    satellite2.pos = satellite2.pos + (satellite2.speed * dt)
    satellite3.pos = satellite3.pos + (satellite3.speed * dt)




# Step #2 - Representation of the Earth and the Moon using classes

from __future__ import division, print_function
from vpython import *

scene=canvas(title="Moon & Earth")

class Planet: #we define the new class for planets
    def __init__ (self, mass, radius, position, color): #assign the input parameters of the planet
        self.mass = mass
        self.radius = radius
        self.position = position
        self.s = sphere (pos=position, radius=radius, color=color)

# We call the classes and create the earth and the moon representations
earth = Planet (5.98e24,6.4e6,vector(0,0,0),color.blue)
moon = Planet (0.007346e24,1.736e6,vector(0,378e6+6.4e6,0),color.white)




# Step #3 - Representation of three satellites (LEO,MEO,GEO) orbiting around Earth, using planet and satellite classes.
# It also contains the enhanced version, where the satellites have wings and cones that represent the surface area covered of the earth.

from __future__ import division, print_function
from vpython import *

scene=canvas(title="Satellite")

class Planet: #we define the new class for planets
    def __init__ (self, mass, radius, position, color): #assign the input parameters of the planet
        self.mass = mass
        self.radius = radius
        self.position = position
        self.s = sphere (pos=position, radius=radius, color=color,texture={'file':textures.earth})

class Satellite: #we define the new class for planets
    def __init__ (self, planet, speed, position, size, color, coverage_angle): #assign the input parameters of the satellite
        self.planet = planet.position
        self.speed = speed
        self.position = position
        self.size = size
        self.color = color
        self.angle = coverage_angle # The coverage angle indicates the angle of the cone that covers the surface of the Earth
        self.s = sphere (pos=self.position, radius=self.size, color=self.color, make_trail=True, retain=80)
        
        # Representation of the wings of the satellite and the cone
        self.wing = box (pos=self.position, axis=vector(self.position.x,0,0), length=self.size*8, height=self.size/5, width=self.size*2/3, up=vector(0,0,1))
        self.cone = cone(pos=self.planet, axis=self.position, radius=(tan(self.angle*pi/360)*(position.y)),color=self.color,opacity=0.6)

    def updatePosition(self,dt):
        r = self.planet - self.position
        F = (G * earth_mass * norm(r)) / (mag2(r)) # F/Unidad de Masa
        self.speed = self.speed + F*dt
        self.position = self.position + self.speed*dt
        self.s.pos = self.position
        self.wing.pos = self.position
        self.cone.axis = self.position
       
G = 6.67e-11
earth_position = vector (0,0,0)
earth_mass = 5.98e24
earth = Planet (earth_mass,6.4e6,earth_position,color.white)

# We define first the LEO satellite    
LEO_speed = vector((sqrt(G*earth_mass/(6.4e6+163e3))),0,0)
LEO_position = vector(0,(6.4e6+163e4),0)
LEO = Satellite (earth,LEO_speed,LEO_position,1e6,color.red,60) # We call the class satellite to plot the 3D Representation

# We define then the MEO satellite 
MEO_speed = vector((sqrt(G*earth_mass/(6.4e6+20e6))),0,0)
MEO_position = vector(0,(6.4e6+20e6),0)
MEO = Satellite (earth,MEO_speed,MEO_position,1e6,color.yellow,12) # We call the class satellite to plot the 3D Representation

# Finally we define the GEO satellite 
GEO_speed = vector((sqrt(G*earth_mass/(6.4e6+36e6))),0,0)
GEO_position = vector(0,(6.4e6+36e6),0)
GEO = Satellite (earth,GEO_speed,GEO_position,1e6,color.green,17) # We call the class satellite to plot the 3D Representation

dt=10
while True:
    rate(100)
    LEO.updatePosition (dt)
    MEO.updatePosition (dt)
    GEO.updatePosition (dt)





# Step #4 - The program represents a constellation of 4 orbits with 3 satellites each, 
# all orbiting the Earth in orbits crossing over the poles

from __future__ import division, print_function
from vpython import *

scene=canvas(title="Satellite Constellation")
    

class Planet: #we define the new class for planets
    def __init__ (self, mass, radius, position, color): #assign the input parameters of the planet
        self.mass = mass
        self.radius = radius
        self.position = position
        self.s = sphere (pos=position, radius=radius, color=color,texture={'file':textures.earth})

class Satellite: #we define the new class for planets
    def __init__ (self, planet, speed, position, size, color, coverage_angle): #assign the input parameters of the Satellite
        self.planet = planet.position
        self.speed = speed
        self.position = position
        self.size = size
        self.color = color
        self.angle = coverage_angle
        self.s = sphere (pos=self.position, radius=self.size, color=self.color, make_trail=True, retain=80)
                
        self.wing = box (pos=self.position, axis=vector(self.position.x,0,0), length=self.size*8, height=self.size/5, width=self.size*2/3, up=vector(0,0,1))
        self.cone = cone(pos=self.planet, axis=self.position, radius=(tan(self.angle*pi/360)*(mag(self.planet-self.position))),color=self.color,opacity=0.5)

    def updatePosition(self,dt):
        r = self.planet - self.position
        F = (G * earth_mass * norm(r)) / (mag2(r)) # F/Unidad de Masa
        self.speed = self.speed + F*dt
        self.position = self.position + self.speed*dt
        self.s.pos = self.position
        self.wing.pos = self.position
        self.cone.axis = self.position

G = 6.67e-11
earth_position = vector (0,0,0)
earth_mass = 5.98e24
earth = Planet (earth_mass,6.4e6,earth_position,color.white)

# Definition of the original MEO  
MEO_speed_original = (sqrt(G*earth_mass/(6.4e6+20e6)))*vector(1,0,0)
MEO_position_original = vector(0,(6.4e6+20e6),0)

# From the original MEO, we will define the rotation vectors
number_orbit = 0
number_satellite = 0
rotation_orbit = [0,pi/4,pi/2,3*pi/4] # This vector rotates along "y" axis, to account for the 4 different orbits over the poles
rotation_satellite = [0,120*pi/180,240*pi/180] # This vector rotates along "z" axis, to account for the 3 different satllites per orbit
MEO = [0] # We create an empty list to save all the created satellites

#  Now we create the loop that goes through every rotation angle 
#  (first all the orbits, then all satellites in each orbit)
for angle_orbit in rotation_orbit:
    number_orbit += 1 
    #print ("Orbit: ",number_orbit)
    colour = [0,color.red,color.yellow,color.blue,color.green]
    for angle_satellite in rotation_satellite:
        number_satellite += 1
        #print ("Sat: ",number_satellite)
        # We define the axis around which the satellites will rotate (they're perpendicular to the orbit plane)
        rotation_axis = [0,vector(0,0,1),vector(-0.5,0,-0.5),vector(-1,0,0),vector(0.5,0,-0.5)]
        # We first calculate the rotated position
        MEO_position = rotate (MEO_position_original,angle=angle_orbit,axis=vector(0,1,0))
        MEO_position = rotate (MEO_position,angle=angle_satellite,axis=rotation_axis[number_orbit])
        # And then we compute the rotated speed
        MEO_speed = rotate (MEO_speed_original,angle=angle_orbit,axis=vector(0,1,0))
        MEO_speed = rotate (MEO_speed,angle=angle_satellite,axis=rotation_axis[number_orbit])
        # Finally, we compute the satellite in 3D with the class previously deifned
        MEO.append (Satellite (earth,MEO_speed,MEO_position,1e6,colour[number_orbit],12))
        
# Loop to update de position of all the satellites
dt=10
while True:
    rate(100)
    MEO[1].updatePosition (dt)
    MEO[2].updatePosition (dt)
    MEO[3].updatePosition (dt)
    MEO[4].updatePosition (dt)
    MEO[5].updatePosition (dt)
    MEO[6].updatePosition (dt)
    MEO[7].updatePosition (dt)
    MEO[8].updatePosition (dt)
    MEO[9].updatePosition (dt)
    MEO[10].updatePosition (dt)
    MEO[11].updatePosition (dt)
    MEO[12].updatePosition (dt)
    
    
    

