class Node:

    def __init__(self, x1, x2, y1, y2):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.xCentre = x2-x1
        self.yCentre = y2-y1
        self.env = 0
        self.children = []
        self.ID = 0
   
    # provided node becomes child
    def set_child(self, child):
        self.children.append(child)
    
    def set_env(self, env):
        self.env = env

    # returns the environment of the node
    # if the node has no environment it searchs the right kid for coordinates and calls itself again
    def get_env(self, x, y):
        if self.env == 0:
            square = self.get_square(x, y)
            return self.children[square].get_env(x, y)
        return self.env

    # used by get_env(), returns the square of the coordinates
    def get_square(self, x, y):
        square = 0
        if x >= self.xCentre:
            square |= 2
        if y >= self.yCentre:
            square |= 1
        return square