class Node:

    def __init__(self, x1, x2, y1, y2, z1, z2):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.z1 = z1
        self.z2 = z2
        self.xCentre = x2-x1
        self.yCentre = y2-y1
        self.zCentre = z2-z1
        self.env = 0
        self.children = []
        self.ID = 0

    #returns the outer variabels of the node
    def get_outer_var(self):
        return {"x1": self.x1, "x2": self.x2, "y1": self.y1, "y2": self.y2, "z1":  self.z1, "z2": self.z2}
    
    #provided node becomes child
    def set_child(self, child):
        self.children.append(child)
    
    def set_env(self, env):
        self.env = env



    def set_unique_number(self, number):
        self.ID = number

    #returns the environment of the node
    #if the node has no environment it searchs the right kid for coordinates and calls itself again
    def get_env(self, x, y, z):
        if self.env == 0:
            octant = self.get_octant(x, y, z)
            return self.children[octant].get_env(x, y, z)
        return self.env

    def get_octant(self, x, y, z):
        octant = 0
        if x >= self.xCentre:
            octant |= 4
        if y >= self.yCentre:
            octant |= 2
        if z >= self.zCentre:
            octant |= 1
        return octant