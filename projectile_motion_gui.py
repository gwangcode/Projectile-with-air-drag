import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLineEdit, QPushButton, QTextEdit
from PyQt5.QtCore import Qt

class ProjectileWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Projectile Motion Simulator")
        self.setGeometry(100, 100, 1200, 800)

        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        input_widget = QWidget()
        input_layout = QFormLayout()
        self.v0_input = QLineEdit("20")
        self.theta_input = QLineEdit("45")
        self.h0_input = QLineEdit("0")
        self.hf_input = QLineEdit("0")
        self.k_input = QLineEdit("0.01")  # Reduced drag constant
        self.n_input = QLineEdit("2")     # Changed to quadratic drag
        input_layout.addRow("Initial speed v0 (m/s):", self.v0_input)
        input_layout.addRow("Initial angle theta (degrees):", self.theta_input)
        input_layout.addRow("Initial height h0 (m):", self.h0_input)
        input_layout.addRow("Final height hf (m):", self.hf_input)
        input_layout.addRow("Air drag constant k:", self.k_input)
        input_layout.addRow("Exponent n:", self.n_input)
        self.simulate_button = QPushButton("Simulate")
        self.simulate_button.clicked.connect(self.run_simulation)
        input_layout.addWidget(self.simulate_button)
        input_widget.setLayout(input_layout)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        input_layout.addRow("Results:", self.output_text)

        self.figure, self.axes = plt.subplots(3, 1, figsize=(8, 12))
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(input_widget)
        main_layout.addWidget(self.canvas)

    def simulate_motion(self, v0, theta, h0, hf, k, n, include_drag=True):
        g = 9.81
        dt = 0.0001  # Smaller time step for better accuracy
        x = 0.0
        y = h0
        vx = v0 * np.cos(theta)
        vy = v0 * np.sin(theta)
        t = 0.0

        xs = [x]
        ys = [y]
        ts = [t]
        vxs = [vx]
        vys = [vy]

        while y >= hf:  # Changed to >= to handle exact landing
            v = np.sqrt(vx**2 + vy**2)
            if include_drag:
                if v < 1e-10:
                    ax = 0.0
                    ay = -g
                else:
                    drag_factor = k * v**(n - 1)  # Adjusted for correct drag force units
                    ax = -drag_factor * vx
                    ay = -g - drag_factor * vy
            else:
                ax = 0.0
                ay = -g

            x += vx * dt
            y += vy * dt
            vx += ax * dt
            vy += ay * dt
            t += dt

            xs.append(x)
            ys.append(y)
            ts.append(t)
            vxs.append(vx)
            vys.append(vy)

        # Interpolate to find exact point where y = hf
        if y < hf:
            dy = ys[-1] - ys[-2]
            if abs(dy) > 1e-10:
                frac = (hf - ys[-2]) / dy
                x_final = xs[-2] + frac * (xs[-1] - xs[-2])
                t_final = ts[-2] + frac * (ts[-1] - ts[-2])
                vx_final = vxs[-2] + frac * (vxs[-1] - vxs[-2])
                vy_final = vys[-2] + frac * (vys[-1] - vys[-2])
                y_final = hf
            else:
                x_final = xs[-1]
                t_final = ts[-1]
                vx_final = vxs[-1]
                vy_final = vys[-1]
                y_final = ys[-1]
        else:
            x_final = xs[-1]
            t_final = ts[-1]
            vx_final = vxs[-1]
            vy_final = vys[-1]
            y_final = ys[-1]

        xs[-1] = x_final
        ys[-1] = y_final
        ts[-1] = t_final
        vxs[-1] = vx_final
        vys[-1] = vy_final

        return np.array(xs), np.array(ys), np.array(ts), vxs, vys, x_final, y_final, t_final, vx_final, vy_final

    def run_simulation(self):
        try:
            v0 = float(self.v0_input.text())
            theta_deg = float(self.theta_input.text())
            h0 = float(self.h0_input.text())
            hf = float(self.hf_input.text())
            k = float(self.k_input.text())
            n = float(self.n_input.text())
            theta = np.radians(theta_deg)

            xs, ys, ts, vxs, vys, x_final, y_final, t_final, vx_final, vy_final = self.simulate_motion(v0, theta, h0, hf, k, n, True)
            xs_no_drag, ys_no_drag, ts_no_drag, _, _, x_final_no_drag, y_final_no_drag, t_final_no_drag, vx_final_no_drag, vy_final_no_drag = self.simulate_motion(v0, theta, h0, hf, 0, n, False)

            idx_max = np.argmax(ys)
            t_max = ts[idx_max]
            x_max = xs[idx_max]
            y_max = ys[idx_max]
            idx_max_no_drag = np.argmax(ys_no_drag)
            t_max_no_drag = ts_no_drag[idx_max_no_drag]
            x_max_no_drag = xs_no_drag[idx_max_no_drag]
            y_max_no_drag = ys_no_drag[idx_max_no_drag]

            vf = np.sqrt(vx_final**2 + vy_final**2)
            theta_f = np.degrees(np.arctan2(vy_final, vx_final))
            vf_no_drag = np.sqrt(vx_final_no_drag**2 + vy_final_no_drag**2)
            theta_f_no_drag = np.degrees(np.arctan2(vy_final_no_drag, vx_final_no_drag))

            output = (
                "With Air Drag:\n"
                f"Total time of travel: {t_final:.2f} s\n"
                f"Time at highest point: {t_max:.2f} s\n"
                f"Coordinates of highest point: ({x_max:.2f}, {y_max:.2f})\n"
                f"Final speed: {vf:.2f} m/s\n"
                f"Final angle: {theta_f:.2f} degrees\n"
                f"Coordinates of final point: ({x_final:.2f}, {y_final:.2f})\n\n"
                "Without Air Drag:\n"
                f"Total time of travel: {t_final_no_drag:.2f} s\n"
                f"Time at highest point: {t_max_no_drag:.2f} s\n"
                f"Coordinates of highest point: ({x_max_no_drag:.2f}, {y_max_no_drag:.2f})\n"
                f"Final speed: {vf_no_drag:.2f} m/s\n"
                f"Final angle: {theta_f_no_drag:.2f} degrees\n"
                f"Coordinates of final point: ({x_final_no_drag:.2f}, {y_final_no_drag:.2f})"
            )
            self.output_text.setText(output)

            self.axes[0].clear()
            self.axes[1].clear()
            self.axes[2].clear()
            self.axes[0].plot(xs, ys, label='With Air Drag')
            self.axes[0].plot(xs_no_drag, ys_no_drag, '--', label='No Air Drag')
            self.axes[0].set_xlabel('x (m)')
            self.axes[0].set_ylabel('y (m)')
            self.axes[0].set_title('y vs x')
            self.axes[0].grid(True)
            self.axes[0].legend()

            self.axes[1].plot(ts, ys, label='With Air Drag')
            self.axes[1].plot(ts_no_drag, ys_no_drag, '--', label='No Air Drag')
            self.axes[1].set_xlabel('t (s)')
            self.axes[1].set_ylabel('y (m)')
            self.axes[1].set_title('y vs t')
            self.axes[1].grid(True)
            self.axes[1].legend()

            self.axes[2].plot(ts, xs, label='With Air Drag')
            self.axes[2].plot(ts_no_drag, xs_no_drag, '--', label='No Air Drag')
            self.axes[2].set_xlabel('t (s)')
            self.axes[2].set_ylabel('x (m)')
            self.axes[2].set_title('x vs t')
            self.axes[2].grid(True)
            self.axes[2].legend()

            self.figure.tight_layout()
            self.canvas.draw()

        except ValueError:
            self.output_text.setText("Error: Please enter valid numerical values.")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = ProjectileWindow()
    window.show()
    sys.exit(app.exec_())
