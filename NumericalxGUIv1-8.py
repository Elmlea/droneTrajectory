import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def projectile_trajectory(v0, m, initial_y, dt, c, A):
    # Constants
    p = 1.45661  # Air density (kg/m^3)
    k = c * p * A / 2

    g = 9.81  # Acceleration due to gravity (m/s^2)

    # Initial conditions
    x, y = 0.0, initial_y
    vx, vy = v0, 0.0

    # Lists to store the trajectory
    x_values, y_values = [x], [y]

    # Euler method for numerical integration
    while y >= 0:
        v = np.sqrt(vx**2 + vy**2)
        ax = -k / m * v * vx
        ay = -g - k / m * v * vy

        x += vx * dt
        y += vy * dt

        vx += ax * dt
        vy += ay * dt

        x_values.append(x)
        y_values.append(y)

    # Calculate kinetic energy at y = 0
    kinetic_energy = 0.5 * m * vy**2

    return x_values, y_values, kinetic_energy

def plot_trajectory(x_values, y_values):
    fig, ax = plt.subplots()
    ax.plot(x_values, y_values)
    ax.set_xlabel('Horizontal Distance (m)')
    ax.set_ylabel('Vertical Height (m)')
    ax.set_title('Projectile Trajectory')
    return fig

def calculate_trajectory():
    # Get input from the GUI
    v0 = float(entry_v0.get())
    m = float(entry_mass.get())
    initial_y = float(entry_initial_y.get())
    c = float(entry_drag_coefficient.get())
    A = float(entry_cross_sectional_area.get())
    dt = 0.001  # Time step for numerical integration (s)

    # Calculate projectile trajectory until y = 0
    x_values, y_values, kinetic_energy = projectile_trajectory(v0, m, initial_y, dt, c, A)

    # Display results in the GUI
    result_label.config(text=f"Final x position: {x_values[-1]:.2f} m\n"
                             f"Final y position: {y_values[-1]:.2f} m\n"
                             f"Horizontal velocity (v_x) at y = 0: {x_values[-1]/dt:.2f} m/s\n"
                             f"Vertical velocity (v_y) at y = 0: {y_values[-1]/dt:.2f} m/s\n"
                             f"Kinetic energy at y = 0: {kinetic_energy:.2f} J")

    # Create the plot
    fig = plot_trajectory(x_values, y_values)

    # Update the plot in the GUI
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

# Create the GUI
root = Tk()
root.title("Projectile Trajectory Calculator")

# Input fields
label_v0 = Label(root, text="Initial Velocity (m/s):")
label_v0.pack()
entry_v0 = Entry(root)
entry_v0.pack()

label_mass = Label(root, text="Mass (kg):")
label_mass.pack()
entry_mass = Entry(root)
entry_mass.pack()

label_initial_y = Label(root, text="Initial Vertical Position (m):")
label_initial_y.pack()
entry_initial_y = Entry(root)
entry_initial_y.pack()

label_drag_coefficient = Label(root, text="Drag Coefficient:")
label_drag_coefficient.pack()
entry_drag_coefficient = Entry(root)
entry_drag_coefficient.pack()

label_cross_sectional_area = Label(root, text="Cross-sectional Area (m^2):")
label_cross_sectional_area.pack()
entry_cross_sectional_area = Entry(root)
entry_cross_sectional_area.pack()

# Calculate button
calculate_button = Button(root, text="Calculate", command=calculate_trajectory)
calculate_button.pack()

# Results
result_label = Label(root, text="", justify=LEFT)
result_label.pack()

# Plot frame
plot_frame = Frame(root)
plot_frame.pack()

root.mainloop()
