# Interactive polygon plot
# Simple version of iplotpoly - Includes subdivide and rotate buttons

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox

# Interactive polygon plot
def iplotpoly(inputA):
    global vX, vY, plot_center, text_message
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

    # Rotate polyline 90 degrees clockwise
    def rotateClockwiseCallback(event):
        global vX, vY

        tempX = np.array(vX)
        tempY = np.array(vY)

        tempX = tempX - plot_center[0];
        tempY = tempY - plot_center[1];

        tempX2 = tempY;
        tempY2 = -tempX;

        tempX2 = tempX2 + plot_center[0];
        tempY2 = tempY2 + plot_center[1];

        vX = tempX2.tolist();
        vY = tempY2.tolist();

        redrawPlot(vX, vY)

        return


    # Subdivide polyline
    def subdivideCallback(event):
        global vX, vY

        if (len(vX) == 0): return

        tempX = [ vX[0] ]
        tempY = [ vY[0] ]

        i = 0;
        while (i+1 < len(vX)):
            tempX.append((vX[i] + vX[i+1])/2)
            tempY.append((vY[i] + vY[i+1])/2)

            tempX.append(vX[i+1])
            tempY.append(vY[i+1])

            i = i+1

        vX = tempX;
        vY = tempY;
        nump = len(vX)

        redrawPlot(vX, vY)

        return



    # Enable button
    def enable_button(button, buttonCallback):

        enabled_color = 'white'
        enabled_hover_color = 'green'

        button.on_clicked(buttonCallback)
        button.color = enabled_color
        button.hovercolor = enabled_hover_color

        return


    # Enable buttons by linking to callbacks
    def enable_all_buttons():
        global subdivide_button
        global rotateCW_button

        enable_button(subdivide_button, subdivideCallback)
        enable_button(rotateCW_button, rotateClockwiseCallback)

        # BUG: Button color does not change until the mouse is moved.
        # Other programmers have complained about this bug/feature.

        plt.draw()

        return


    # After polyline vertices are created or read, draw polyline
    #   and enable buttons and pick event
    def drawPolylineEnableButtonsPickEvent(vX, vY):
        global cid_pick

        plt.subplot(111)

        outputPlotMessage("Select and move vertices")
        drawPolyline(vX,vY)
        plt.draw()

        cid_pick = fig.canvas.mpl_connect('pick_event', pickPointCallback)

        # Enable buttons
        enable_all_buttons()

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

            drawPolylineEnableButtonsPickEvent(vX, vY)

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


    # Create buttons
    def createButtons(left_pos, button_width, button_height):
        global ax_rotateCW_button, rotateCW_button
        global ax_subdivide_button, subdivide_button
        global ax_save_button, save_button

        disabled_color = 'dimgray'

        ax_subdivide_button = plt.axes([left_pos, 0.8, button_width, button_height])
        subdivide_button = Button(ax_subdivide_button, "Subdivide")
        subdivide_button.color = disabled_color
        subdivide_button.hovercolor = disabled_color

        ax_rotateCW_button = plt.axes([left_pos, 0.7, button_width, button_height])
        rotateCW_button = Button(ax_rotateCW_button, "Rotate CW")
        rotateCW_button.color = disabled_color
        rotateCW_button.hovercolor = disabled_color

        return



    # Initialize global variables
    text_message = None
    iv_selected = 0
    plot_center = [ 0.5, 0.5 ]
    default_output_filename = 'controlpts.txt'

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

    createButtons(0.91, 0.08, 0.05)

    # User interactively selects vertices
    cid_addVertexCallback = fig.canvas.mpl_connect('button_press_event', addVertexCallback)

    plt.show()

    return
