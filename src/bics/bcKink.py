'''
    Module to define the boundary condition of the problem:
        - u(t, 0)   = m
        - u(t, Nxt) = -m 
'''
def bcI(ic0=0, ic1=0, V=0):
    '''
        B.C. for the initial point
    '''
    import numpy as np 
    import inc
    #----------------------
    m = inc.vp2
    return -m
#---------------------------------------
def bcF(ic0=0, ic1=0, V=0):
    '''
        B.C. for the final point
    '''
    import numpy as np 
    import inc
    #----------------------
    m = inc.vp2
    return m
