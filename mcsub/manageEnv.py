from environment import environment
import math
class manageEnv:

    def __init__(self):
        self.envArray = []      #all environments that are not default
        self.envDefault =  environment("define in manageEnv", 1, 1, 1, 1, "x+y+z=0", 0, 1)   #default environment

    def add_environment(self, name, mua, mus, excitAnisotropy, refractiveIndex, height, radius, default, bins):
        if default == 0:
            height = [height[0]*bins, height[1]*bins]
            radius *= bins

            self.envArray.append(environment(name, mua, mus, excitAnisotropy, refractiveIndex, height, radius, default))
        else:
            if self.envDefault.name != "define in manageEnv":
                raise ValueError("Multiple default environments defined")
            self.envDefault = environment(name, mua, mus, excitAnisotropy, refractiveIndex, height, radius, default)
    

    #find the corresponding environement for a point
    #brute-force method, which is very inefficent
    #currently not used in this program
    def find_env(self, radius, height):
        for env in self.envArray:
            if env.in_environment(radius, height):
                return env.get_variables()
        return self.envDefault.get_variables()
    
    #updates a 1D environment lista
    #based on the radius and height only spefic bins are checked if they are in the corresponding environment
    def assign_env(self, envList, pixel):
        for env in self.envArray:
            for h in range(math.ceil(env.height[0]), math.ceil(env.height[1])+1):
                for r in range(0, math.ceil(env.radius)+1):
                    if env.in_environment(r+0.5, h+0.5):
                        envList[int(abs(h)*pixel+ abs(r))] = env.get_variables()
        return envList
    
    #returns the variables of the default environment
    def get_default_variables(self):
        return self.envDefault.get_variables()