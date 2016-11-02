# IMPORT STUFF
import profile
import pstats
import os
# END IMPORT


os.system('rm tmpprofile')
profile.run('from moments import moments; moments()','tmpprofile')

p = pstats.Stats('tmpprofile')
p.sort_stats('cumulative').print_stats(10)

