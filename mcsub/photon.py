import math

class photon:

    def __init__(self, nr):
        self.id = nr
        self.x = 0
        self.y = 0
        self.z = 0
        self.weight = 1.0
        self.status = 1

    #update the postional variables of the photon
    def update_positon(self, ux, uy, uz):
        self.x += ux
        self.y += uy
        self.z += uz
    
    def get_position(self):
        return [self.x, self.y, self.z]
    
    def get_spherical_position(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
    
    def get_cylindrical_position(self):
        return math.sqrt(self.x*self.x + self.y*self.y)
    
    def get_planar_position(self):
        return abs(self.z)
    
    #update weight of the photon
    def update_weight(self, w):
        self.weight = w

    def get_weight(self):
        return self.weight

    #photon dies
    def update_status(self): 
        self.status = 0
    
    def get_status(self):
        return self.status