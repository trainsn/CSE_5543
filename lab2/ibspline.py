# Interactive polygonal line plot.
# Interactively set polyline vertices or read vertices from a file.
# Buttons to add control points
# Button to save vertices (control points) to a file.

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Interactive polygon plot
def iplotpoly(mode, inputA):
    global vX, vY, plot_center, text_message, save_message
    global jv, iv_selected, numv
    global default_output_filename, filename_textbox
    global cid_bpress, cid_pick
    global save_message_xloc
    global ndegree
    global knots
    global add_status

    # Draw a polygonal line through points vX to vY
    def drawPolyline(vX, vY):
        plt.plot(vX, vY, color="magenta")
        plt.plot(vX, vY, 'ob', picker=True, pickradius=5)
        global numv, ndegree
        global knots
        k = ndegree + 1

        for j in range(k-1, numv):
            q = np.zeros((k, 50, 2))
            u = np.linspace(0, 1, 50) * (knots[j+1] - knots[j]) + knots[j]
            for i in range(j-k+1, j+1):
                q[i - (j-k+1), :, 0] = vX[i]
                q[i - (j-k+1), :, 1] = vY[i]
            for r in range(1, k):
                h = k - r
                # construct control points for order h spline
                for i in range(j, j-h, -1):
                    alpha = (u - knots[i]) / (knots[i+h] - knots[i])
                    q[i - (j-k+1)] = (1.0 - alpha)[:, np.newaxis] * q[i-1 - (j-k+1)] + alpha[:, np.newaxis] * q[i - (j-k+1)]
            plt.plot(q[k-1, :, 0], q[k-1, :, 1], 'b-')

        return

    # Add a control point by Boehm algorithm
    def addControlPoint(iv):
        global vX, vY, numv, ndegree
        global knots
        if iv < ndegree:
            outputPlotMessage("Select and move vertices")
            outputSaveMessage("For the selected control point q_h, h should be at least d")
            return
        q = np.zeros((numv + 1, 2))
        k = ndegree + 1
        # Keep first j ??? k + 1 control points
        q[:iv-k+2, 0] = np.array(vX[:iv-k+2])
        q[:iv-k+2, 1] = np.array(vY[:iv-k+2])
        # Keep last n ??? j + 1 control points
        q[iv+1:, 0] = np.array(vX[iv:])
        q[iv+1:, 1] = np.array(vY[iv:])
        # Replace k ??? 2 control points from q_{j???k+2} to q_{j???1} with k ??? 1 new control points
        t = (knots[iv] + knots[iv+1]) / 2.
        i = np.arange(iv-k+2, iv+1, 1)
        alpha = t - np.array(knots)[i] / (np.array(knots)[i+k-1] - np.array(knots)[i])
        q[i, 0] = (1 - alpha) * np.array(vX)[i-1] + alpha * np.array(vX)[i]
        q[i, 1] = (1 - alpha) * np.array(vY)[i-1] + alpha * np.array(vY)[i]
        vX, vY = q[:, 0].tolist(), q[:, 1].tolist()
        numv += 1
        knots.insert(iv + 1, t)
        redrawPlot(vX, vY)

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
        plt.xlim([xL, xR])
        plt.ylim([yB, yT])

        drawPolyline(vX, vY)
        outputPlotMessage("Select and move vertices")
        clearSaveMessage()
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

    # Set a message related to saving control points to file
    def outputSaveMessage(s):
        global save_message, save_message_xloc
        text_yloc=1.02
        if save_message is not None:
            save_message.remove()

        save_message=plt.text(save_message_xloc-0.14, text_yloc, s, transform=ax.transAxes)
        return

    # Set the save message to ''
    def clearSaveMessage():
        global save_message

        if save_message is not None:
            save_message.remove()
            save_message = None
        return

    # CallbackS

    # add control points
    def addcpCallback(event):
        global vX, vY, numv
        global add_status

        outputPlotMessage("Select a existing control point for control pointing adding task")
        clearSaveMessage()
        add_status = True

        return

    # Create TextBox containing output filename and
    #   output message text box
    def createSaveTextBoxes(default_filename):
        global filename_textbox, save_message_xloc

        ax_filename_textbox = fig.add_axes([save_message_xloc, 0.94, 0.3, 0.03])
        filename_textbox = TextBox(ax_filename_textbox, 'Output filename:', default_filename)
        plt.draw()

        return

    # Save control points to a file
    def saveCallback(event):
        global vX, vY
        global knots

        if (len(vX) == 0): return

        filename = filename_textbox.text;

        if (filename == '' or filename.isspace()):
            # Should pop up an error window/message, but the GUI code would become
            #   even more complicated.
            outputSaveMessage('Illegal filename. Reset filename to default. No file saved.')
            filename_textbox.set_val(default_output_filename)
            plt.draw()

        else:
            # Save file
            f = open(filename, "w")
            f.write("BSPLINE\n{:d} {:d} {:d}\n".format(2, numv, ndegree))
            for knot in knots:
                f.write("{:.2f} ".format(knot))
            f.write("\n")
            for i in range(numv):
                f.write("{:f} {:f}\n".format(vX[i], vY[i]))

            outputSaveMessage('Control points saved to file ' + filename)
            plt.draw()

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
        global addcp_button
        global rotateCW_button
        global save_button
        global default_output_filename

        enable_button(addcp_button, addcpCallback)
        enable_button(save_button, saveCallback)

        # BUG: Button color does not change until the mouse is moved.
        # Other programmers have complained about this bug/feature.

        createSaveTextBoxes(default_output_filename)

        plt.draw()

        return

    # After polyline vertices are created or read, draw polyline
    #   and enable buttons and pick event
    def drawPolylineEnableButtonsPickEvent(vX, vY):
        global cid_pick

        plt.subplot(111)

        # Store current axes limits
        xL, xR = ax.get_xlim()
        yB, yT = ax.get_ylim()

        # Bug: Clearing the axes resets the navigation toolbar so you can't return
        #  to previous states.
        # In particular, you can't return to the starting (0,0), (1,1) configuration.
        plt.cla()

        # Restore axes limits
        plt.xlim([xL, xR])
        plt.ylim([yB, yT])

        outputPlotMessage("Select and move vertices")
        drawPolyline(vX, vY)
        plt.draw()

        cid_pick = fig.canvas.mpl_connect('pick_event', pickPointCallback)

        # Enable buttons
        enable_all_buttons()

        return

    # Add a vertex at (event.xdata, event.ydata)
    def addVertexCallback(event):
        global jv, vX, vY, cid_bpress

        plt.subplot(111)

        # Store current axes limits
        xL, xR = ax.get_xlim()
        yB, yT = ax.get_ylim()

        # Restore axes limits
        plt.xlim([xL, xR])
        plt.ylim([yB, yT])

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
        global add_status

        if len(event.ind) < 1:
            return

        if (event.ind[0] >= 0) and (event.ind[0] < len(vX)):
            iv_selected = event.ind[0]
            if not add_status:
                cid_moveCallback = fig.canvas.mpl_connect('button_release_event', moveVertexCallback)
            else:
                addControlPoint(iv_selected)
                add_status = False
        return

    # Create buttons
    def createButtons(left_pos, button_width, button_height):
        global ax_rotateCW_button, rotateCW_button
        global ax_addcp_button, addcp_button
        global ax_save_button, save_button

        disabled_color = 'dimgray'

        ax_addcp_button = plt.axes([left_pos, 0.8, button_width, button_height])
        addcp_button = Button(ax_addcp_button, "Add CP")
        addcp_button.color = disabled_color
        addcp_button.hovercolor = disabled_color

        ax_rotateCW_button = plt.axes([left_pos, 0.7, button_width, button_height])
        rotateCW_button = Button(ax_rotateCW_button, "Rotate CW")
        rotateCW_button.color = disabled_color
        rotateCW_button.hovercolor = disabled_color

        ax_save_button = plt.axes([left_pos, 0.6, button_width, button_height])
        save_button = Button(ax_save_button, "Save")
        save_button.color = disabled_color
        save_button.hovercolor = disabled_color

        return

    # load vX and vY from a textfile
    def loadTextFile(filename):
        global vX, vY, jv, numv, ndegree
        global knots

        vX, vY = [], []
        with open(filename, "r") as fh:
            assert(next(fh).strip("\r\n") == "BSPLINE")
            for line in fh:
                if line[0] != "#":
                    ndim, numv, ndegree = [int(tmp) for tmp in line.split()]
                    break
            assert(ndim == 2)
            knots = [float(tmp) for tmp in next(fh).split()]
            for line in fh:
                x, y = [float(tmp) for tmp in line.split()]
                vX.append(x)
                vY.append(y)

        assert(numv == len(vX))
        jv = numv

        return

    # Initialize global variables
    text_message = None
    save_message = None
    filename_textbox = None
    save_message_xloc = 0.5
    iv_selected = 0
    plot_center = [ 0.5, 0.5 ]
    default_output_filename = 'controlpts.txt'

    # If inputA is an integer, set number of vertices to inputA
    # If inputA is a string, read control points from inputA
    if mode == "integer":
        numv, ndegree = [int(item) for item in inputA.split(",")]
        jv = 0
        vX = []
        vY = []
        k = ndegree + 1
        knots = [0. for i in range(k)]
        knots.extend([float(i-k+1) for i in range(k, numv)])
        knots.extend([float(numv-k+1) for i in range(numv, numv+k)])
    elif mode == "file":
        loadTextFile(inputA)
    else:
        print('Only support mode integer and file')
        print('Exiting.')
        return

    if (numv < 2):
        print('Illegal number of vertices: ', numv)
        print('Requires at least two vertices.')
        print('Exiting')
        return

    add_status = False

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])

    plt.subplots_adjust(right=0.9)

    outputSelectVertexLocation(0)

    # get plot lower left (LL) and upper right (UR)
    axLL = ax.transData.transform((0, 0))
    axUR = ax.transData.transform((1, 1))

    createButtons(0.91, 0.08, 0.05)

    if (jv < numv):
        # User interactively selects vertices
        cid_addVertexCallback = fig.canvas.mpl_connect('button_press_event', addVertexCallback)
    else:
        # Vertices already read in from input file
        drawPolylineEnableButtonsPickEvent(vX, vY)

    plt.show()

    return

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("--mode", required=True, type=str,
                        help="program mode: integer or file")
parser.add_argument("--inputA", required=True, type=str,
                        help="input to the file")

if __name__ == '__main__':
    args = parser.parse_args()
    iplotpoly(args.mode, args.inputA)
