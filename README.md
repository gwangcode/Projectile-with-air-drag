Projectile Motion Simulator
This repository contains a Python application that simulates projectile motion with and without air drag, visualized in a graphical user interface (GUI) using PyQt5 and Matplotlib. The simulator allows users to input initial conditions and parameters, plots the motion trajectories, and displays key results such as flight time, highest point coordinates, and final conditions.
Features

Input initial speed, angle, heights, air drag constant, and exponent.
Simulate projectile motion with customizable air drag (f = k * v^n).
Display three plots: y vs x, y vs t, and x vs t, comparing trajectories with and without air drag.
Output total flight time, time at the highest point, coordinates of the highest point, final speed, final angle, and final coordinates.
GUI interface built with PyQt5 for user interaction.

Requirements

Python 3.x
Libraries:
numpy
matplotlib
PyQt5


PyInstaller (for creating the executable)

Installation

Clone the repository:git clone https://github.com/gwangcode/Projectile-with-air-drag
cd projectile-motion-simulator


Create a virtual environment (optional but recommended):python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate


Install dependencies:pip install -r requirements.txt

(Note: Create a requirements.txt file with the following content if not already present:numpy
matplotlib
PyQt5
pyinstaller



Usage

Run the script directly:python projectile_motion_gui_fixed.py


Enter the initial conditions in the GUI:
Initial speed (v0 in m/s)
Initial angle (theta in degrees)
Initial height (h0 in m)
Final height (hf in m)
Air drag constant (k)
Exponent (n)


Click "Simulate" to generate plots and results.
View the plots and output in the GUI window.

Building the Executable
To create a standalone .exe file for Windows:

Ensure PyInstaller is installed (pip install pyinstaller).
Run the following command in the terminal:pyinstaller --onefile --windowed --hidden-import=matplotlib --hidden-import=PyQt5 projectile_motion_gui.py


Find the executable in the dist folder.

Contributing
Feel free to fork this repository, submit issues, or send pull requests for improvements or bug fixes.
License
This project is licensed under the MIT License - see the LICENSE file for details.
Acknowledgments

Built with help from the xAI community and Grok AI assistance.
Utilizes open-source libraries: NumPy, Matplotlib, and PyQt5.
