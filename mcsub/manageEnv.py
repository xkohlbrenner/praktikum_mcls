from environment import environment

class manageEnv:

    def __init__(self):
        self.envArray = []      #all environments that are not default
        self.envDefault =  environment("define in manageEnv", 1, 1, 1, "x+y+z=0", 1)   #default environment

    def addEnvironment(self, name, mua, mus, excitAnisotropy, formula, default):
        if default == 0:
            self.envArray.append(environment(name, mua, mus, excitAnisotropy, formula, default))
        else:
            if self.envDefault.name != "define in manageEnv":
                raise ValueError("Multiple default environments defined")
            self.envDefault = environment(name, mua, mus, excitAnisotropy, formula, default)
    
    def findEnv(self, ix, iy, iz):
        for env in self.envArray:
            if env.in_environment(ix, iy, iz):
                return env.get_variables()
        return self.envDefault.get_variables()