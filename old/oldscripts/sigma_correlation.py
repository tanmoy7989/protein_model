#!/usr/bin/env python
import os, sys, numpy as np
import Go as Go, utils

PdbName = sys.argv[1]
BBType = sys.argv[2]
PolyResName = sys.argv[3]

# template for script to run relative entropy optimization
srel_str = '''
#!/usr/bin/env python
import os, sys, numpy as np
import Go as Go

Go.InPdb = %(NATIVEPDB)s
AATraj = %(AATRAJ)s
Go.Prefix = %(PREFIX)s
Go.Seq = ['%(POLYRESNAME)s'] * 15 
Go.Temps = np.logspace(np.log10(%(TLOW)g, %(THIGH)g, 8)
Go.TempSet = Go.Temps[np.argmin( abs(Go.Temps - %(TREF)g) )]


