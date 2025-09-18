import tkinter as tk
from tkinter import ttk, messagebox
import subprocess
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os
import sys

curpath = 'C:\\Users\\chehp\\OneDrive\Desktop\\all\\Git-Projects\\Stiff_equation\\'
# fnames = "C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/pFile.txt", "C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/ansFile.txt", "C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/gTolFile.txt", "C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/rTolFile.txt"
fnames = curpath+'data\\pFile.txt', curpath+'data\\ansFile.txt', curpath+'data\\gTolFile.txt', curpath+'data\\rTolFile.txt'
tnames = "RK", "Answer", "Global", "Local"

class TrajectoryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Numerical Trajectories Viewer")
        self.root.geometry("900x700")

        style = ttk.Style()
        style.configure("TLabel", font=("Segoe UI", 11))
        style.configure("TEntry", font=("Segoe UI", 11))
        style.configure("TButton", font=("Segoe UI", 11), padding=6)
        style.configure("TCheckbutton", font=("Segoe UI", 11))

        param_frame = ttk.LabelFrame(root, text="Simulation Parameters", padding=10)
        param_frame.pack(side="top", fill="x", padx=10, pady=10)

        self.entries = {}
        params = [
            ("Step", "0.005"),
            ("Error", "0"),
            ("X Left", "-1000"),
            ("X Right", "1000"),
            ("Y Lower", "-1000"),
            ("Y Upper", "1000"),
            ("T Left", "0"),
            ("T Right", "0.1"),
        ]

        for i, (param, default) in enumerate(params):
            ttk.Label(param_frame, text=param).grid(row=i // 2, column=(i % 2) * 2, sticky="e", padx=5, pady=5)
            entry = ttk.Entry(param_frame, width=12)
            entry.insert(0, default)
            entry.grid(row=i // 2, column=(i % 2) * 2 + 1, padx=5, pady=5)
            self.entries[param] = entry


        btn_frame = ttk.Frame(param_frame)
        btn_frame.grid(row=4, column=0, columnspan=4, pady=10)

        ttk.Button(btn_frame, text="Run Simulation", command=self.run_simulation).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Exit", command=sys.exit).pack(side="left", padx=5)


        traj_frame = ttk.LabelFrame(root, text="Trajectories", padding=10)
        traj_frame.pack(side="left", fill="y", padx=10, pady=5)

        self.trajectory_vars = []
        for i in range(4):
            var = tk.IntVar(value=1 if i == 0 else 0) 
            cb = ttk.Checkbutton(traj_frame, text=tnames[i], variable=var, command=self.plot_trajectories)
            cb.pack(anchor="w", pady=2)
            self.trajectory_vars.append(var)

        plot_frame = ttk.Frame(root)
        plot_frame.pack(side="right", expand=True, fill="both", padx=10, pady=10)

        self.fig, self.ax = plt.subplots(figsize=(7, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.canvas.mpl_connect("scroll_event", self.on_scroll)

        self.data = {}

    def run_simulation(self):
        try:
            step = float(self.entries["Step"].get())
            error = float(self.entries["Error"].get())
            x_left = float(self.entries["X Left"].get())
            x_right = float(self.entries["X Right"].get())
            y_lower = float(self.entries["Y Lower"].get())
            y_upper = float(self.entries["Y Upper"].get())
            t_left = float(self.entries["T Left"].get())
            t_right = float(self.entries["T Right"].get())
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numbers.")
            return

        try:
            subprocess.run([
                curpath+'x64\\Release\\KSR1.exe ',
                "0",
                "7",
                "13",
                f"{step}",
                f"{error}",
                f"{t_left}",
                f"{x_left}",
                f"{y_lower}",
                f"{t_right}",
                f"{x_right}",
                f"{y_upper}"
            ], check=True)
        except Exception as e:
            messagebox.showerror("Execution Error", str(e))
            return

        self.data.clear()
        for i in range(4):
            fname = fnames[i]
            if os.path.exists(fname):
                arr = np.loadtxt(fname)
                if arr.ndim == 1: 
                    arr = arr.reshape(1, -1)
                self.data[i] = arr
        self.plot_trajectories()

    def plot_trajectories(self):
        self.ax.clear()
        colors = [(1., 0., 0.), (0., 1., 0.), (0., 0., 1.), (1., 0., 1.)]
        for i, var in enumerate(self.trajectory_vars):
            if var.get() and i in self.data:
                arr = self.data[i]
                color = colors[i]
                self.ax.plot(arr[:, 0], arr[:, 1], color=color, linewidth = 1, label=tnames[i]+', x1')
                self.ax.scatter(arr[:, 0], arr[:, 1], color=color, s = 25)
                newcolor = color[0] / 2, color[1] / 2, color[2] / 2
                self.ax.plot(arr[:, 0], arr[:, 2], color=newcolor, linewidth = 1, label=tnames[i]+', x2')
                self.ax.scatter(arr[:, 0], arr[:, 2], color=newcolor, s = 25)
        self.ax.set_xlabel("t")
        self.ax.set_ylabel("x1, x2")
        self.ax.legend()
        self.ax.grid(True, linestyle="--", alpha=0.5)
        self.canvas.draw()

    def on_scroll(self, event):
        base_scale = 1.2
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        xdata = event.xdata 
        ydata = event.ydata 

        if xdata is None or ydata is None:
            return

        if event.button == "up":
            scale_factor = 1 / base_scale
        elif event.button == "down":
            scale_factor = base_scale
        else:
            scale_factor = 1

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

        relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])

        self.ax.set_xlim([xdata - (1 - relx) * new_width, xdata + relx * new_width])
        self.ax.set_ylim([ydata - (1 - rely) * new_height, ydata + rely * new_height])
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = TrajectoryApp(root)
    root.mainloop()