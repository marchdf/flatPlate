#
# Post-process data using VisIt. Generates data files of
# - the wall shear stress
# - the pressure force in x, y, z
#
# Usage:
#    visit -nowin -cli -s /path/to/pp_aux_vars.py
#

# ========================================================================
#
# Imports
#
# ========================================================================
import sys
import os

# ========================================================================
#
# Function definitions
#
# ========================================================================


def save_curve(plotnum, fname):
    """Save curve data"""
    SetActivePlots(plotnum)
    HideActivePlots()
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = "."
    SaveWindowAtts.fileName = fname
    SaveWindowAtts.family = 0
    SaveWindowAtts.format = SaveWindowAtts.CURVE
    SaveWindowAtts.width = 1024
    SaveWindowAtts.height = 1024
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.None
    SaveWindowAtts.forceMerge = 0
    # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
    HideActivePlots()


# ========================================================================
#
# Main
#
# ========================================================================


# Open files
fdir = os.path.abspath('results')
fname = os.path.join(fdir, 'flatPlate.e')
OpenDatabase("localhost:" + fname, 0)

# Wall shear stress
AddPlot("Curve", "operators/Lineout/tau_wall", 1, 1)
LineoutAtts = LineoutAttributes()
LineoutAtts.point1 = (0, -0.5, 0)
LineoutAtts.point2 = (2, -0.5, 0)
LineoutAtts.interactive = 0
LineoutAtts.ignoreGlobal = 0
LineoutAtts.samplingOn = 0
LineoutAtts.numberOfSamplePoints = 100
LineoutAtts.reflineLabels = 0
SetOperatorOptions(LineoutAtts, 1)

# Pressure force x
AddPlot("Curve", "operators/Lineout/pressure", 1, 1)
LineoutAtts = LineoutAttributes()
LineoutAtts.point1 = (0, -0.5, 0)
LineoutAtts.point2 = (2, -0.5, 0)
LineoutAtts.interactive = 0
LineoutAtts.ignoreGlobal = 0
LineoutAtts.samplingOn = 0
LineoutAtts.numberOfSamplePoints = 100
LineoutAtts.reflineLabels = 0
SetOperatorOptions(LineoutAtts, 1)

DrawPlots()

# Hide all curve plots
SetActiveWindow(1)
SetActivePlots(0)
HideActivePlots()
SetActivePlots(1)
HideActivePlots()

# Go to final time
nfinal = TimeSliderGetNStates() - 1
SetTimeSliderState(nfinal)

# Save the wall shear stress
save_curve(0, os.path.join(fdir, "tau_wall"))

# Save the pressure forces
save_curve(1, os.path.join(fdir, "pressure"))

sys.exit()
