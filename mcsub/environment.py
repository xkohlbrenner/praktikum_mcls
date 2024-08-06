class environment:

    def __init__(self, name, mua, mus, excitAnisotropy, refractiveIndex, height, radius, default):
        self.name = name
        self.mua = mua
        self.mus = mus
        self.excitAnisotropy = excitAnisotropy
        self.refIndex = refractiveIndex
        self.default = default                      #1, if default env, else 0
        if self.default == 0:
            self.height = height
            self.radius = radius

    
    # boolean function returns true, if coordinates are in the environment
    def in_environment(self, radius, height):
        if self.default == 1:
            return True
        return (radius <= self.radius and self.height[0] <= height <= self.height[1])

    # return the variables of the environment
    def get_variables(self):
        return {"name": self.name, "mua": self.mua, "mus": self.mus, "excitAnisotropy": self.excitAnisotropy, "refIndex": self.refIndex}