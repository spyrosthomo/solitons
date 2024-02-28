'''
  
'''
class Potential:
    import inc 
    #----------------------------
    def __init__(self, phin, l, m, null1=0):
        import numpy as np 
        #-------
        self.potentialName = "phi4"
        self.potentialForm = f"V = l(m-phi^2)^2"
        self.forFileNames  = f"phi4_l{l}_m{m}"
        self.saveLoc       = "../figures/.pdf"
        self.xti           = Potential.inc.xti
        self.xtf           = Potential.inc.xtf
        self.Nxt           = Potential.inc.Nxt
        self.l             = l 
        self.m             = m
        self.phin          = phin
        self.xt            = np.linspace(self.xti, self.xtf, self.Nxt+1)
    #-----------------------------
    def potential(self):
        import numpy as np 
        m    = self.m
        l    = self.l
        phin = self.phin
        #V    = l*(m-phin**2)**2
        Vp   = 4*l*(m-phin**2)*phin
        return Vp 
    #----------------------------

