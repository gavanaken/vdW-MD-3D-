from Tkinter import *
import sys
from MD3D import *
import os


currentSystem = None
root = Tk()
root.resizable(0,0)
windowbg = 'slate gray'
def disable_event():
    pass
root.protocol("WM_DELETE_WINDOW", disable_event)
root.configure(background=windowbg)
root.title("VDW simulator")
# Define a function to close the window.
def quit(event=None):
    root.destroy()
    sys.exit()

# Cause pressing <Esc> to close the window.
root.bind('<Escape>', quit)



# -----------------------------------------------------------------------------
# Create a frame within the main window.
# -----------------------------------------------------------------------------
# The following few frame objects will contain the widgets needed to run a simulation.
frame = Frame(root, background=windowbg)
frame.pack()
# Create a "terminal-style" textbox to offer updates
termVar = StringVar(frame)
T = Label(frame, textvariable=termVar, bg='black', fg='white', height=4, width=60)
T.config(borderwidth=4)
T.pack(side='left')
termVar.set('TERMINAL\nset your configuration parameters below...')

def getError():
    print('got here')
    if currentSystem.dtError == True:
        print('got here too')
        #global term
        T.configure(text='TERMINAL\nlooks like your particles moved too fast\nlower your dt')

e = Button(frame, text='get Error', command=getError, height=4)
e.pack(side='left')

# Define an input variable and add an entry box so the user can change its value.
namePrompt = Label(frame, text="            name this run: ", background=windowbg)
namePrompt.pack(side='left')

runName = StringVar()
runName.set('')
name_input = Entry(frame, width=30, textvariable=runName)
name_input.pack(side='left')


systemSetup = Label(root, text="\nset parameters for this run:", font="Helvetica 12 bold", background=windowbg)
systemSetup.pack(side='top')
setup = Frame(root, background=windowbg)
setup.pack(side='top')


# Input number of particles:
npartVar = StringVar(setup)
particleChoices = [2, 8, 27, 64, 125, 216, 343, 512, 729, 1000]
npartVar.set(216)

partPrompt = Label(setup, text="     # of particles:", background=windowbg)
partPrompt.pack(side='left')
npartDrop = OptionMenu(setup, npartVar, *particleChoices)
npartDrop.config(highlightthickness=0)
npartDrop.pack(side='left')

# Input boxlength:
boxlengthVar = StringVar(setup)
boxlengthChoices = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
boxlengthVar.set(4)

boxlengthPrompt = Label(setup, text="     boxlength (per side):", background=windowbg)
boxlengthPrompt.pack(side='left')
boxlengthDrop = OptionMenu(setup, boxlengthVar, *boxlengthChoices)
boxlengthDrop.pack(side='left')
boxlengthDrop.config(highlightthickness=0)

# Input temperature:
tempVar = StringVar(setup)
tempChoices = [1, 1.095, 1.5, 5, 10, 20, 50, 100, 150, 200, 300]
tempVar.set(1)

tempPrompt = Label(setup, text="     temperature:", background=windowbg)
tempPrompt.pack(side='left')
tempDrop = OptionMenu(setup, tempVar, *tempChoices)
tempDrop.pack(side='left')
tempDrop.config(highlightthickness=0)


# Input rcut:
rcutVar = StringVar(setup)
rcutChoices = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
rcutVar.set(2.5)

rcutPrompt = Label(setup, text="     radial cutoff:", background=windowbg)
rcutPrompt.pack(side='left')
rcutDrop = OptionMenu(setup, rcutVar, *rcutChoices)
rcutDrop.pack(side='left')
rcutDrop.config(highlightthickness=0)

# Input dt:
dtVar = StringVar(setup)
dtChoices = [0.01, 0.001, 0.0005, 0.0001, 0.00008, 0.00006, 0.00004, 0.00002, 0.000015, 0.00001, 0.000001]
dtVar.set(0.0001)

dtPrompt = Label(setup, text="     timestep(dt):", background=windowbg)
dtPrompt.pack(side='left')
dtDrop = OptionMenu(setup, dtVar, *dtChoices)
dtDrop.pack(side='left')
dtDrop.config(highlightthickness=0)

designSetup = Label(root, text="\n\nset design for this run:", font="Helvetica 12 bold", background=windowbg)
designSetup.pack(side='top')
design = Frame(root, background=windowbg)
design.pack(side='top')

# Particle size:
sizeVar = StringVar(setup)
sizeChoices = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300]
sizeVar.set(60)

sizePrompt = Label(design, text="     particle size:", background=windowbg)
sizePrompt.pack(side='left')
sizeDrop = OptionMenu(design, sizeVar, *sizeChoices)
sizeDrop.pack(side='left')
sizeDrop.config(highlightthickness=0)

# Coloring option:
copVar = StringVar(setup)
copChoices = ['standard', 'force', 'velocity']
copVar.set('standard')

copPrompt = Label(design, text="     coloring option:", background=windowbg)
copPrompt.pack(side='left')
copDrop = OptionMenu(design, copVar, *copChoices)
copDrop.pack(side='left')
copDrop.config(highlightthickness=0)


# Particle colors:
colorVar = StringVar(setup)
colorChoices = ['inferno', 'binary', 'gray', 'bone', 'Blues', 'Greys', 'Oranges', 'Greens', 'RdPu', 'PuBuGn', 'spring',
                'autumn', 'winter', 'summer', 'hot', 'random']
colorVar.set('random')

colorPrompt = Label(design, text="     particle colors:", background=windowbg)
colorPrompt.pack(side='left')
colorDrop = OptionMenu(design, colorVar, *colorChoices)
colorDrop.pack(side='left')
colorDrop.config(highlightthickness=0)


# Background color:
bgVar = StringVar(setup)
bgChoices = ['black', 'white', 'whitesmoke', 'antiquewhite', 'gray', 'darkgray', 'dimgray', 'slategray',
             'lightslategrey', 'rosybrown', 'lightsteelblue', 'random']
bgVar.set('random')

bgPrompt = Label(design, text="     background color:", background=windowbg)
bgPrompt.pack(side='left')
bgDrop = OptionMenu(design, bgVar, *bgChoices)
bgDrop.pack(side='left')
bgDrop.config(highlightthickness=0)

# Data collection test
dataSetup = Label(root, text="\n\nset data collection parameters:", font="Helvetica 12 bold", background=windowbg)
dataSetup.pack(side='top')
data = Frame(root, background=windowbg)
data.pack(side='top')

# nEquil
nEquil_input = Label(data, text="     equilibration steps: ", background=windowbg)
nEquil_input.pack(side='left')

nEquilVar = StringVar(data)
nEquil_entry = Entry(data, width=8, textvariable=nEquilVar)
nEquil_entry.pack(side='left')
nEquilVar.set(500)

# nRun
nRun_input = Label(data, text="     data collection steps: ", background=windowbg)
nRun_input.pack(side='left')

nRunVar = StringVar(data)
nRun_entry = Entry(data, width=8, textvariable=nRunVar)
nRun_entry.pack(side='left')
nRunVar.set(5000)

dataCol = Frame(root, background=windowbg)
dataCol.pack(side='top')
filler = Label(dataCol, width=1, text="", background=windowbg)
filler.pack(side="top")

# nData
nData_input = Label(dataCol, text="collect thermo-data every", background=windowbg)
nData_input.pack(side='left')

nDataVar = StringVar(dataCol)
nData_entry = Entry(dataCol, width=8, textvariable=nDataVar)
nData_entry.pack(side='left')
nDataVar.set(50)

nData_input2 = Label(dataCol, text="steps", background=windowbg)
nData_input2.pack(side='left')

# nCoords
nCoord_input = Label(dataCol, text="                             capture coordinate positions every", background=windowbg)
nCoord_input.pack(side='left')

nCoordVar = StringVar(dataCol)
nCoord_entry = Entry(dataCol, width=8, textvariable=nCoordVar)
nCoord_entry.pack(side='left')
nCoordVar.set(100)

nData_input2 = Label(dataCol, text="steps", background=windowbg)
nData_input2.pack(side='left')

def executeAni():
    if colorVar.get() == 'random':
        cmap = cm.get_cmap(random.choice(colorChoices[:-1]))
    else:
        cmap = cm.get_cmap(colorVar.get())
    if bgVar.get() == 'random':
        bgcolor = random.choice(bgChoices[:-1])
    else:
        bgcolor = bgVar.get()
    colors = cmap(numpy.linspace(0, 1, npartVar.get()))
    cop = copVar.get()
    global currentSystem
    currentSystem = genSystem(int(npartVar.get()), float(boxlengthVar.get()),
                 float(tempVar.get()), float(rcutVar.get()), float(dtVar.get()), cmap, int(sizeVar.get()), bgcolor, cop)
    setTerm("TERMINAL\nrunning animation...")
    runAnimation(currentSystem)

def stopAni():
    global currentSystem
    setTerm('TERMINAL\nset your configuration parameters below...')
    currentSystem.setExit()

def setTerm(message):
    termVar.set(message)

def _collectData():
    PElist, KElist, TElist, coordlist, Plist, Glist = runSimulation(currentSystem, currentSystem.dt, int(nRunVar.get()),
                                                      int(nDataVar.get()), int(nCoordVar.get()))
    setTerm("TERMINAL\nSuccess! files saved to .../data as" + "\n" + runName.get() + "_PE/_KE/_TE/_coord")

    direct = os.path.dirname(__file__)

    PEfile = open(direct+"/data/"+name_input.get()+"_PE", "w")
    PEfile.write(str(PElist))

    KEfile = open(direct+"/data/"+name_input.get()+'_KE', 'w')
    KEfile.write(str(KElist))

    TEfile = open(direct+"/data/"+name_input.get() + '_TE', 'w')
    TEfile.write(str(TElist))

    coordfile = open(direct+"/data/"+name_input.get() + '_coord', 'w')
    coordfile.write(str(coordlist))

    Pfile = open(direct + "/data/" + name_input.get() + '_Pressure', 'w')
    Pfile.write(str(Plist))

    Gfile = open(direct + "/data/" + name_input.get() + '_GInf', 'w')
    Gfile.write(str(Glist))


def collectData():
    setTerm("TERMINAL\nbeginning run and generating data...")
    return _collectData()

def _runEquil():
    global currentSystem
    currentSystem = genSystem(int(npartVar.get()), float(boxlengthVar.get()),
                              float(tempVar.get()), float(rcutVar.get()), float(dtVar.get()), None, int(sizeVar.get()),
                              None, None)
    currentSystem = equilibrate(currentSystem, int(nEquilVar.get()))
    setTerm("TERMINAL\nequilibration complete\nready to collect data")
    getData.config(text= ' collect data ', command=collectData)

def runEquil():
    setTerm("TERMINAL\nequilibrating...")
    return _runEquil()

def serial_test():
    press = []
    TE = []
    global currentSystem
    direct = os.path.dirname(__file__)
    serialfile = open(direct + "/data/" + name_input.get() + '_serial', 'w')
    lengths = [12, 10, 8, 6, 4]
    for size in lengths:
        currentSystem = genSystem(125, size, 1, 2.5, 0.001, None, None, None, None)
        currentSystem = equilibrate(currentSystem, 500)
        PElist, KElist, TElist, coordlist, Plist, Glist = runSimulation(currentSystem, 0.001, 500, 10, 1000)
        nodtTE = []
        for element in TElist:
            nodtTE.append(element[1])
        TE.append([size,numpy.mean(nodtTE)])
        nodtP = []
        for element in Plist:
            nodtP.append(element[1])
        press.append([size,numpy.mean(nodtP)])
        serialfile.write(str(size)+" TE: " + str(TElist))
        serialfile.write(str(size)+" P: " + str(Plist))
    serialfile.write('alltogether: ' + " TE: " + str(TE))
    serialfile.write('alltogether: ' + " P: " + str(press))



buttons = Frame(root, background=windowbg)
buttons.pack(side='bottom')
filler = Label(buttons, width=1, text="", background=windowbg)
filler.pack(side="top")
# Create a button to perform the calculation and pack it into the frame.
stopAnimate = Button(buttons, text=' stop animation ', command=stopAni, padx=4, bg='firebrick1')
stopAnimate.pack(side='left')
filler = Label(buttons, width=1, text="", background=windowbg)
filler.pack(side="left")
doAnimate = Button(buttons, text=' animate ', command=executeAni, padx=4, bg='green2')
doAnimate.pack(side='left')
filler = Label(buttons, width=1, text="", background=windowbg)
filler.pack(side="left")
getData = Button(buttons, text=' equilibrate ', command=runEquil, padx=4, bg='DeepSkyBlue2')
getData.pack(side='left')
filler = Label(buttons, width=1, text="", background=windowbg)
filler.pack(side="left")
leaveButt = Button(buttons, text="exit", command=quit, padx=10, bg = 'gold')
leaveButt.pack(side='left')
#testButt = Button(buttons, text="testSerial", command=serial_test)
#testButt.pack(side="left")

# -----------------------------------------------------------------------------
# Activate the window.
# -----------------------------------------------------------------------------
root.mainloop()