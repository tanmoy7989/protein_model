
#!/usr/bin/env python

''' simple class template to supply configs to the model
    '''

import numpy as np
from const import *

nice = lambda x : 0.0 if x is None else x

# create a config object and populate with defaults
class Config(object):
    def __init__(self):
        for k,v in DEFAULTS.iteritems():
            self.__setattr__(k, v)
       
    def __repr__(self):
        s1 = HELPSTR
        s2 = 'Current config:\n'
        for k, v in DEFAULTS.iteritems(): s2 += '%s: %s\n' % (k, self.__dict__[k])
        return s1 + '\n' + s2

    

