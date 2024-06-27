import sympy as sp
class environment:

    def __init__(self, name, mua, mus, excitAnisotropy, height, radius, default):
        self.name = name
        self.mua = mua
        self.mus = mus
        self.excitAnisotropy = excitAnisotropy
        self.default = default                      #1, if default env, else 0
        if self.default == 0:
            self.height = height
            self.radius = radius

    

    def in_environment(self, radius, height):
        if self.default == 1:
            return True
        return (radius <= self.radius and self.height[0] <= height <= self.height[1])

    
    def get_variables(self):
        return {"name": self.name, "mua": self.mua, "mus": self.mus, "excitAnisotropy": self.excitAnisotropy}