from environment import environment

class manageEnv:

    def __init__(self):
        self.envArray = []      #all environments that are not default
        self.envDefault =  environment("define in manageEnv", 1, 1, 1, "x+y+z=0", 0, 1)   #default environment

    def addEnvironment(self, name, mua, mus, excitAnisotropy, formula, space, default):
        if default == 0:
            self.envArray.append(environment(name, mua, mus, excitAnisotropy, formula, space, default))
        else:
            if self.envDefault.name != "define in manageEnv":
                raise ValueError("Multiple default environments defined")
            self.envDefault = environment(name, mua, mus, excitAnisotropy, formula, space, default)
    
    def find_env(self, ix, iy, iz):
        for env in self.envArray:
            if env.in_environment(ix, iy, iz):
                return env.get_variables()
        return self.envDefault.get_variables()
    
    def assign_env(self, envList, pixel):
        for env in self.envArray:
            for ix in range(env.space["x"][0], env.space["x"][1]+1):
                for iy in range(env.space["y"][0], env.space["y"][1]+1):
                    for iz in range(env.space["z"][0], env.space["z"][1]+1):
                        if env.in_environment(ix+0.5, iy+0.5, iz+0.5):
                            envList[int(abs(ix)*pixel*pixel+abs(iy)*pixel+abs(iz)*pixel)] = env.get_variables()
        return envList
    
    def get_default_variables(self):
        return self.envDefault.get_variables()