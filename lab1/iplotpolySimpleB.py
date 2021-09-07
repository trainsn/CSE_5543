# Interactive polygon plot
# Simple version of iplotpoly - Move vertices after creating polyline.

import sys
import numpy as np
import matplotlib.pyplot as plt

# Interactive polygon plot
def iplotpoly(inputA):
    global vX, vY, text_message
    global jv, iv_selected, numv
    global cid_bpress, cid_pick


    # Draw a polygonal line through points vX to vY
    def drawPolyline(vX, vY):
        plt.plot(vX, vY, color="magenta")
        plt.plot(vX, vY, 'ob', picker=True, pickradius=5)
        return


    # Redraw the polyline plot
    def redrawPlot(vX, vY):
        plt.subplot(111)

        # Store current axes limits
        xL, xR = ax.get_xlim()
        yB, yT = ax.get_ylim()

        # Bug: Clearing the axes resets the navigation toolbar so you can't return
        #  to previous states.
        # In particular, you can't return to the starting (0,0), (1,1) configuration.
        plt.cla()

        # Restore axes limits
        plt.xlim([xL,xR])
        plt.ylim([yB,yT])

        drawPolyline(vX, vY)
        outputPlotMessage("Select and move vertices")
        plt.draw()

        return


    # Output a message above the plot
    def outputPlotMessage(s):
        global text_message
        text_yloc=1.02
        if text_message is not None:
            text_message.remove()

        text_message=plt.text(0.0, text_yloc, s, transform=ax.transAxes)
        return


    # Output the message: Select vertex str(j) location
    def outputSelectVertexLocation(j):
        s = "Select vertex " + str(j) + " location"
        outputPlotMessage(s)
        return


    # CallbackS


    # After polyline vertices are created or read, draw polyline
    #   and enable pick event
    def drawPolylineEnablePickEvent(vX, vY):
        global cid_pick

        plt.subplot(111)

        outputPlotMessage("Select and move vertices")
        drawPolyline(vX,vY)
        plt.draw()

        cid_pick = fig.canvas.mpl_connect('pick_event', pickPointCallback)

        return


    # Add a vertex at (event.xdata, event.ydata)
    def addVertexCallback(event):
        global jv, vX, vY, cid_bpress

        plt.subplot(111)
        x=event.xdata
        y=event.ydata

        # Check that x and y are defined
        if (x is None) or (y is None): return

        # check that event.x and event.y are within draw region
        if ((event.x < axLL[0]) or (event.x > axUR[0])): return
        if ((event.y < axLL[1]) or (event.y > axUR[1])): return

        if (jv < numv):
            vX.append(x)
            vY.append(y)
            plt.plot(x, y, 'ob')
            jv = jv+1

        if (jv < numv):
            outputSelectVertexLocation(jv)
            plt.draw()
        else:
            # Disconnect button press callback before creating pick event
            fig.canvas.mpl_disconnect(cid_addVertexCallback)

            drawPolylineEnablePickEvent(vX, vY)

        return


    # Move selected vertex to location (event.xdata, event.ydata)
    def moveVertexCallback(event):
        global cid_moveCallback, iv_selected, vX, vY
        global points, polyline

        plt.subplot(111)
        x=event.xdata
        y=event.ydata

        # Check that x and y are defined
        if (x is None) or (y is None): return

        # check that event.x and event.y are within draw region
        if ((event.x < axLL[0]) or (event.x > axUR[0])): return
        if ((event.y < axLL[1]) or (event.y > axUR[1])): return

        vX[iv_selected] = x
        vY[iv_selected] = y
        redrawPlot(vX, vY)

        fig.canvas.mpl_disconnect(cid_moveCallback)


    # Pick and move a vertex
    def pickPointCallback(event):
        global iv_selected, cid_moveCallback

        if (len(event.ind) < 1): return

        if (event.ind[0] >= 0) and (event.ind[0] < len(vX)):
            iv_selected = event.ind[0]
            cid_moveCallback = fig.canvas.mpl_connect('button_release_event', moveVertexCallback)

        return


    # Initialize global variables
    text_message = None
    iv_selected = 0

    # If inputA is an integer, set number of vertices to inputA
    # If inputA is a string, read control points from inputA
    if (isinstance(inputA, int)):
        numv = inputA
        jv = 0
        vX = []
        vY = []
    else:
        print('Illegal input: ', inputA)
        print('Exiting.')
        return

    if (numv < 2):
        print('Illegal number of vertices: ', numv)
        print('Requires at least two vertices.')
        print('Exiting')
        return

    fig, ax = plt.subplots(1,1,figsize=(12,8))
    plt.xlim([0.0,1.0])
    plt.ylim([0.0,1.0])

    plt.subplots_adjust(right=0.9)

    outputSelectVertexLocation(0)

    # get plot lower left (LL) and upper right (UR)
    axLL = ax.transData.transform((0,0));
    axUR = ax.transData.transform((1,1));

    # User interactively selects vertices
    cid_addVertexCallback = fig.canvas.mpl_connect('button_press_event', addVertexCallback)

    plt.show()

    return
