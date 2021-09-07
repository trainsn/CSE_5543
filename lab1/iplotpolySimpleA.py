# Interactive polygon plot
# Simple version of iplotpoly - Select vertex locations and draw polyline.

import sys
import numpy as np
import matplotlib.pyplot as plt


# Interactive polygon plot
def iplotpoly(inputA):
    global vX, vY, text_message
    global jv, numv
    global cid_bpress


    # Draw a polygonal line through points vX to vY
    def drawPolyline(vX, vY):
        plt.plot(vX, vY, color="magenta")

        # picker set to True for more sophisticated version of iplotpoly
        #   where vertices can be moved.
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
        outputPlotMessage("Done")
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
            # Disconnect button press callback
            fig.canvas.mpl_disconnect(cid_addVertexCallback)

            redrawPlot(vX, vY)

        return


    # Initialize global variables
    text_message = None
    iv_selected = 0
    plot_center = [ 0.5, 0.5 ]

    # If inputA is an integer, set number of vertices to inputA
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
