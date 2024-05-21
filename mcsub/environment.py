import sympy as sp

class environment:

    def __init__(self, name, mua, mus, excitAnisotropy, formula, default):
        self.name = name
        self.mua = mua
        self.mus = mus
        self.excitAnisotropy = excitAnisotropy
        self.default = default                      #1, if default env, else 0
        if self.default == 0:
            self.formula = sp.sympify(formula)          #formula to describe the environment, input has to be a string

    

    def in_environment(self, ix, iy, iz):
        if self.default == 1:
            return True
        x, y, z = sp.symbols('x y z')
        return self.formula.subs({x: ix, y: iy, z: iz})
    
    def get_variables(self):
        return {"mua": self.mua, "mus": self.mus, "excitAnisotropy": self.excitAnisotropy}