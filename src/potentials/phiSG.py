'''
  
'''
class Potential:
    import inc 
    #----------------------------
    def __init__(self, phin, l=0, m=0, null1=0):
        import numpy as np 
        #-------
        self.potentialName = "phiSH"
        self.potentialForm = f"V = 1-cos(phi)"
        self.forFileNames  = f"phi4SG"
        self.saveLoc       = f"../figures/{self.forFileNames}.pdf"
        self.xti           = Potential.inc.xti
        self.xtf           = Potential.inc.xtf
        self.Nxt           = Potential.inc.Nxt
        self.phin          = phin
        self.xt            = np.linspace(self.xti, self.xtf, self.Nxt+1)
    #-----------------------------
    def potential(self):
        import numpy as np 
        phin = self.phin
        V    = 1-np.cos(phin)
        return V 
    #----------------------------

