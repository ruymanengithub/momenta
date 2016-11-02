#! /usr/bin/env python

from Tkinter import *
from pdb import set_trace as stop
from copy import copy

varnames = ['CLTHRESH','CLMINAREA','clsubradial']

def updateEntry(self):
    print '\n'
    for n in varnames:
        self.values[n] = self.Vars[n].get()
        # self.Vars[n].set(self.values[n])
        print n,self.values[n]
    print '\n'
    return

def closeWindow(self):
    self.win.destroy()
    return

def makeWindow(self):
    
    self.win = Tk()

    frame1 = Frame(self.win)
    frame1.pack()
    
    for i in range(len(varnames)):
        n = varnames[i]
        Label(frame1, text=n).grid(row=i, column=0, sticky=W)
        self.Vars[n] = StringVar()
        # try: self.Vars[n].set(self.values[n])
        # except KeyError : self.Vars[n].set('-99.0')
        self.Vars[n].set('4')
        Entry(frame1, textvariable= self.Vars[n]).grid(row=i, column=1, sticky=W)
    
    
    frame2 = Frame(self.win)       # Row of buttons
    frame2.pack()
    b1 = Button(frame2,text="Update",command=self.updateEntry)
    b1.pack(side=LEFT);
    b2 = Button(frame2,text="Close",command=self.closeWindow)
    b2.pack(side=LEFT);
    
    self.win.mainloop()
    
    return

class GUI(dict):
    
    from CLgui import updateEntry,closeWindow,makeWindow
    
    def __init__(self,data):
        self.values = copy(data)
        self.win = None
        self.Vars = {}

def CLgui(values):
    import sys
    print '\n HAS BUGS \n'
    sys.exit()
    
    thisGUI = GUI(values)
    thisGUI.makeWindow()
    newvalues = copy(thisGUI.values)
    
    return newvalues
