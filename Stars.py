from visual import *
from visual.graph import *

from math import *
from random import *

class Space:
    def __init__(self, lower_left_back_vector, upper_right_front_vector, population, max_neighboring_pop, parent = None):
        self.llbv = lower_left_back_vector
        self.urfv = upper_right_front_vector
        
        self.pop = population
        self.max_neighboring_pop = max_neighboring_pop

        self.parent = parent
        self.subspaces = []
        self.to_handle = []
        if parent == None:
            self.setup_subspaces()

        #print 'Drawing up space at', self.llbv, self.urfv, "with population", len(self.pop)
        #self.draw_init()

        self.total_mass = sum(p.mass for p in self.pop)
        if self.total_mass != 0:
            self.center_of_mass = vector(0.0, 0.0, 0.0)
            for p in self.pop:
                self.center_of_mass += p.mass * p.pos
            self.center_of_mass /= self.total_mass
        else:
            self.center_of_mass = self.center()
        
    def create_subspaces(self):
        average = self.center()
        llbvs = [vector(x, y, z) for x in [self.llbv.x, average.x] for y in [self.llbv.y, average.y] for z in [self.llbv.z, average.z]] ## in these loops, the x-loop is the _outermost_
        urfvs = [vector(x, y, z) for x in [average.x, self.urfv.x] for y in [average.y, self.urfv.y] for z in [average.z, self.urfv.z]]
        self.subspaces = [Space(llbvs[i], urfvs[i], \
                                population = [], \
                                max_neighboring_pop = self.max_neighboring_pop, \
                                parent = self) \
                          for i in range(8)]
    def sort_pop_into_subspaces(self):
        for p in self.pop:
            which_octant = self.octant_sort(p.pos)
            self.subspaces[which_octant].pop.append(p)
    def setup_subspaces(self):
        #print 'Setting up subspace at', self.llbv, self.urfv, "with population", len(self.pop)
        if len(self.pop) > self.max_neighboring_pop:
            self.create_subspaces()
            self.sort_pop_into_subspaces()
            for ss in self.subspaces:
                ss.setup_subspaces()

    def contains(self, position):
        return (self.llbv.x <= position.x < self.urfv.x) and \
           (self.llbv.y <= position.y < self.urfv.y) and \
           (self.llbv.z <= position.z < self.urfv.z)
    def octant_sort(self, position):
        average = (self.llbv+self.urfv)/2
        return 4*(0 if position.x < average.x else 1) + \
               2*(0 if position.y < average.y else 1) + \
               1*(0 if position.z < average.z else 1)

    def add_to_subspaces(self, p):
        which_octant = self.octant_sort(p.pos)
        subspace = self.subspaces[which_octant]
        subspace.pop.append(p)
        if len(subspace.pop) > subspace.max_neighboring_pop:
            if subspace.subspaces == []:
                subspace.setup_subspaces() ## this automatically adds `p` to all descendents
            else:
                subspace.add_to_subspaces(p)
    def manage(self):
        ## deepest-first recursion to fill `to_handle` (to "listen to complaints from `self`'s children")
        self.to_handle = []
        if self.subspaces != []:
            for ss in self.subspaces:
                ss.manage()
        else: ## we are at the very bottom
            for p in self.pop:
                if not self.contains(p.pos):
                    self.to_handle.append(p)

        new_population = self.pop[:]
        for p in self.to_handle: ## note that this loop will not run if `self` didn't have any subspaces, nor will t run if we are at the very top
            if self.contains(p.pos):
                self.add_to_subspaces(p)
            ## otherwise, if `p` escaped from not only the child but also the parent (`self`)
            else:
                ## then complain to `self.parent`
                assert self.parent != None
                self.parent.to_handle.append(p)
                new_population.remove(p)
        self.pop = new_population

        if self.pop == []:
            for ss in self.subspaces:
                del ss
            self.subspaces = []

        self.total_mass = sum(p.mass for p in self.pop)
        if self.total_mass != 0:
            self.center_of_mass = vector(0.0, 0.0, 0.0)
            for p in self.pop:
                self.center_of_mass += p.mass * p.pos
            self.center_of_mass /= self.total_mass
        else:
            self.center_of_mass = self.center()

    def draw_init(self):
        pos = self.center()
        dimensions = (self.urfv-self.llbv)*1.0
        self.color = choice([color.red, color.green, color.blue])
        self.drawing = box(pos = pos, \
                           length = dimensions.z,
                           height = dimensions.y,
                           width = dimensions.z, \
                           color = self.color, \
                           opacity = 0 if self.pop == [] else 1.0/len(self.pop))

    def draw(self):
        #return
        self.drawing.opacity = 0.05#0 if self.pop == [] else 1.0/len(self.pop)
        for ss in self.subspaces:
            ss.draw()

    def center(self):
        return (self.llbv+self.urfv)/2
    def radius(self):
        return (self.urfv-self.llbv)/2

    def step(self, dt, fineness, parent):
        PE = 0.0
        if self.subspaces == []:
            if self.pop == []: return 0.0
            average_momentum = sum(p.momentum for p in self.pop)/sum(p.mass for p in self.pop)
            max_relative_KE = max(mag2(p.momentum-average_momentum)/p.mass for p in self.pop)
            if max_KE > 5E51:
                small_dt = dt/fineness
                pe = 0.0
                for i in range(fineness):
                    pe = 0.0
                    for p in self.pop:
                        pe += p.interact(parent, theta = 2)
                    for p in self.pop:
                        p.step(small_dt)
                PE += pe
            else:
                for p in self.pop:
                    PE += p.interact(parent, theta = 2)
                for p in self.pop:
                    p.step(dt)
        else:
            for ss in self.subspaces:
                PE += ss.step(dt, fineness, parent)
        return PE
           

################################################

def ball_uniform(radius):
    rtrn = radius*(2*vector(random(), random(), random()) - vector(1.0, 1.0, 1.0))
    while mag(rtrn) >= radius:
        rtrn = radius*(2*vector(random(), random(), random()) - vector(1.0, 1.0, 1.0))
    return rtrn

def flatten(v):
    return vector(v.x, v.y, 0)

 

c = 3E8
G = 6.67E-11
num_stars = 50
max_mass = 2E30 * 1E11/num_stars  ## so total mass is about the same as our milky way
expected_mass = max_mass/2
expected_total_mass = expected_mass*num_stars
galaxy_radius = 0.5E5*1E16 ## 100,000 light year diameter

class star:
    def __init__(self, init_pos, init_mass, init_vel, col=(1.0, 1.0, 1.0)):
        self.pos = init_pos
        self.mass = init_mass
        self.vel = init_vel
        self.momentum = self.mass*self.vel / sqrt(1-mag2(self.vel)/c**2)
        #col = (1.0, (self.mass/max_mass)**2.0, (self.mass/max_mass)**0.5)
        col = col
        self.visual_object = sphere(pos = self.pos, color = col, radius = galaxy_radius/100)

        self.force_on_me = vector(0.0, 0.0, 0.0)
    def draw(self):
        self.visual_object.pos = self.pos

    def KE(self): ## non-relativistic for now
        return mag2(self.momentum)/self.mass

    def step(self, dt):
        self.pos += dt*(self.momentum+0.5*dt*self.force_on_me)/self.mass
        self.momentum += dt*self.force_on_me
        self.force_on_me = vector(0.0, 0.0, 0.0)
    def interact_with_star(self, other):
        vec1to2 = other.pos - self.pos
        GmM = G*self.mass*other.mass
        forceon1 = GmM/mag2(vec1to2) * norm(vec1to2)
        self.force_on_me += 0.5*forceon1
        other.force_on_me += -0.5*forceon1
        #self.momentum += dt*forceon1
        PE = -GmM/mag(vec1to2)
        return PE
        
    def interact(self, space, theta = 1.0):
        c = space.center_of_mass
        r = mag(space.radius())
        if c==self.pos or r/mag(c-self.pos) > theta:
            PE = 0.0
            if space.subspaces != []:
                for ss in space.subspaces:
                    PE += self.interact(ss, theta = theta)
            else:
                for star in space.pop:
                    if star is not self:
                        PE += self.interact_with_star(star)
            return PE
        else:
            vec1to2 = c - self.pos
            GmM = G*self.mass*space.total_mass
            forceon1 = GmM/mag2(vec1to2) * norm(vec1to2)
            self.force_on_me += 0.5*forceon1
            for star in space.pop:
                star.force_on_me += -0.5*forceon1
            #self.momentum += dt*forceon1
            PE = -GmM/mag(vec1to2)
            return PE
           

star_positions1 = [ball_uniform(galaxy_radius)-0.0*vector(galaxy_radius, 0, 0) for i in range(num_stars)]
star_velocities1 = [vector(+0E5, 0, 0)+ball_uniform(1E4) for i in range(num_stars)]#[sqrt(G*(mag(rad)/galaxy_radius)**3*expected_total_mass / mag(rad))  *  ball_uniform(1.0) for rad in star_positions]
stars1 = [star(init_pos = star_positions1[i], init_mass = expected_mass, init_vel = star_velocities1[i], col = (1.0, 1.0, 0.0)) for i in range(num_stars)]

star_positions2 = [ball_uniform(galaxy_radius)+0.0*vector(galaxy_radius, 0, 0) for i in range(num_stars)]
star_velocities2 = [vector(-0E5, 0, 0)+ball_uniform(1E4) for i in range(num_stars)]#[sqrt(G*(mag(rad)/galaxy_radius)**3*expected_total_mass / mag(rad))  *  ball_uniform(1.0) for rad in star_positions]
stars2 = [star(init_pos = star_positions2[i], init_mass = expected_mass, init_vel = star_velocities2[i], col = (0.0, 1.0, 1.0)) for i in range(num_stars)]
stars = stars1+stars2
scene.autoscale = False

size = 3
S = Space(lower_left_back_vector = vector(-size*galaxy_radius, -size*galaxy_radius, -size*galaxy_radius),\
          upper_right_front_vector = vector(size*galaxy_radius, size*galaxy_radius, size*galaxy_radius), \
          population = stars, max_neighboring_pop = 50)

scale = arrow(pos = vector(0, 0, 0), axis = vector(galaxy_radius, 0, 0), shaftwidth = galaxy_radius/60)

def momentum_to_vel(mom, mass):
    return mom/mass #/ sqrt(1 + mag2(mom/mass)/c**2)

## about 13.6 fps

make_scene = lambda title, posx, posy: \
             gdisplay(title=title, \
                      width=600, height=200, \
                      x=posx*600, y=posy*200, \
                      background=(0.8,0.8,0.8), foreground=(0,0,0))

scene2 = make_scene('Energy v. Time', 0, 0)
scene3 = make_scene('Star Radial Position Histogram', 0, 1)

energyK = gcurve(gdisplay = scene2, color = color.yellow)
energyP = gcurve(gdisplay = scene2, color = color.blue)
energyT = gcurve(gdisplay = scene2, color = color.red)

positions = ghistogram(gdisplay = scene3, bins=[i*0.01*size*galaxy_radius for i in range(200)], accumulate=1)
my_histogram = [0 for i in range(200)]
mh = ghistogram(gdisplay = scene3, bins=[i*0.01*galaxy_radius for i in range(200)])

dt = 5E13


t = 0.0
count = 0
KE = sum(mag2(star.momentum)/(2*star.mass) for star in stars)
PE = 0
for star in stars:
    PE += star.interact(S, theta = 2)
for star in stars:
    star.step(dt)
CM_velocity = vector(0.0, 0.0, 0.0)
for star in stars:
    CM_velocity += star.momentum/star.mass / len(stars)
for star in stars:
    star.momentum -= CM_velocity*star.mass
while True:
    if count % 1 == 0:
        for star in stars:
            if not S.contains(star.pos):
                #print "collision!"
                p, l, u = star.pos, S.llbv, S.urfv
                if not (l.x <= p.x < u.x):
                    star.momentum.x *= -1#-0.999
                    star.pos.x = l.x if p.x < l.x else u.x*0.999999 ## !!!!!
                if not (l.y <= p.y < u.y):
                    star.momentum.y *= -1#-0.999
                    star.pos.y = l.y if p.y < l.y else u.y*0.999999 ## !!!!!
                if not (l.z <= p.z < u.z):
                    star.momentum.z *= -1#-0.999
                    star.pos.z = l.z if p.z < l.z else u.z*0.999999 ## !!!!!
##            c = 0
##            while not S.contains(star.pos):
##                star.momentum *= -1
##                star.pos += dt*star.momentum
##                #S = Space(S.center()-2*S.radius(), S.center()+2*S.radius(), S.pop, S.max_neighboring_pop, parent = None)
##                print c, star.pos
##                c += 1
                ## enlarge Ss
                #S = Space(S.center()-2*S.radius(), S.center()+2*S.radius(), S.pop, S.max_neighboring_pop, parent = None)
##            if not S.contains(star.pos): ## bounce
##                #star.momentum *= -1
##                #star.pos += dt*star.momentum
##                star.pos = ball_uniform(galaxy_radius)
        S.manage()
        #S.draw()
##    for star1 in stars:
##        star1.pos += dt*momentum_to_vel(star1.momentum, star1.mass)
##        for star2 in stars:
##            if star1 is star2: break
##            vec1to2 = star2.pos - star1.pos
##            GmM = G*star1.mass*star2.mass
##            forceon1 = GmM/mag2(vec1to2) * norm(vec1to2)
##            star1.momentum += dt*forceon1
##            star2.momentum += -dt*forceon1
##            PE += GmM/mag(vec1to2)
    #for star in stars:
    #    PE += star.interact(S, theta = 2)
    #for star in stars:
    #    star.step(dt)
    old_PE = PE
    PE = 0
    old_KE = KE
    for star in stars:
        PE += star.interact(S, theta = 2)
    for star in stars:
        star.step(dt)
    #PE = S.step(dt, fineness = 10, parent = S)
    KE = sum(mag2(star.momentum)/(2*star.mass) for star in stars)
    momentum_scale = sqrt((old_PE+old_KE-PE)/KE)
    for star in stars:
        star.momentum *= momentum_scale
    KE = sum(mag2(star.momentum)/(2*star.mass) for star in stars)
    
    energyK.plot(pos = (t, KE))
    energyP.plot(pos = (t, PE))
    energyT.plot(pos = (t, KE+PE))
    positions.plot(data = [mag(star.pos) for star in stars], accumulate = 1)
    #for star in stars:
    #    i = int(floor(mag(star.pos) / (0.01*galaxy_radius)))
    #    my_histogram[i] += 1
    if count % 1 == 0:
        for star in stars:
            star.draw()
    count += 1
    t += dt
    CM_velocity = vector(0.0, 0.0, 0.0)
    for star in stars:
        CM_velocity += star.momentum/star.mass / len(stars)
    for star in stars:
        star.momentum -= CM_velocity*star.mass
    #rate(100)
    #print count, t
