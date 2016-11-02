"""Functions and variables related to flags for moments.py"""

allflags = {'NOPETRO': 2**0L, 'MANYPETRO':2**1L, 'NEGPETRO':2**2L,\
'BADRADIAL':2**3L,'ROBUST':2**4L,'EDMETHOD':2**5L,'CUTMETHOD':2**6L,\
'NORADIAL':2**8L,'':2**9L,'':2**10L,'':2**11L,'':2**12L,'':2**13L,\
'BLANK':2**14L,\
'NONCHECKEDRADIAL':2**15L}

def isflagon(allflags,flag):
    if allflags & flag == flag : return True
    else : return False

def addflag(allflags,flag):
    """Adds a flag to a flag variable."""
    # IMPORT STUFF
    from pdb import set_trace as stop
    # END IMPORT
    if type(allflags) != type(flag) : stop()
    allflags = allflags | flag    
    return allflags
    
