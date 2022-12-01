# GPR on a Single Ping of Sonar Data using GPflow

# Packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import gpflow
import tensorflow as tf
from copy import deepcopy
from matplotlib.axes import Axes
from matplotlib.cm import coolwarm
from check_shapes import check_shapes
import tensorflow_probability as tfp
from typing import Optional


# Import Ping Data, store in Nx1 numpy arrays
dir_path = "/home/parisi/PhilRepo/gpr_ms_thesis/data/wiggles_single_ping.csv"
df = pd.read_csv(dir_path)
x_data = np.array([df['X']])
x_data = np.transpose(x_data)
x_data = x_data - np.mean(x_data)
y_data = np.array([df['Z']])
y_data = np.transpose(y_data)
y_data = y_data - np.mean(y_data)
# plt.plot(x_data, y_data, 'kx', mew=2)
# plt.show()

# Create other dummy data
x_dummy = np.array(
    [
        [0.865], [0.666], [0.804], [0.771], [0.147], [0.866], [0.007], [0.026],
        [0.171], [0.889], [0.243], [0.028],
    ]
)
y_dummy = np.array(
    [
        [1.57], [3.48], [3.12], [3.91], [3.07], [1.35], [3.80], [3.82], [3.49],
        [1.30], [4.00], [3.82],
    ]
)


# Set Variables to Run in Model
X = x_data
Y = y_data


# Run GPR Model

# Create Kernel
p_lengthscale = gpflow.Parameter(0.059, trainable=True)
p_variance = gpflow.Parameter(1.15, trainable=False)
p_noise_var = gpflow.Parameter(0.1, trainable=False)

k = gpflow.kernels.SquaredExponential(
    variance = p_variance,
    lengthscales = p_lengthscale,
)

# Create Model
model = gpflow.models.GPR(
    (X, Y),
    kernel = k,
    noise_variance = 0.02,
)
gpflow.set_trainable(model.likelihood, False) #ensures noise_variance is constant
# Train Model
#opt = gpflow.optimizers.Scipy()
#opt.minimize(model.training_loss, model.trainable_variables)

# Predict / Inference (can't predict if nothing is trainable!)
Xnew = np.array([[0.5]])
model.predict_f(Xnew)
model.predict_y(Xnew)


# Confidence Bounds
Xplot = np.linspace(np.round((x_data[0, 0]*1.05),1), np.round(x_data[x_data.size-1,0]*1.05,1), x_data.size*3)[:, None]
f_mean, f_var = model.predict_f(Xplot, full_cov=False)
y_mean, y_var = model.predict_y(Xplot)
f_lower = f_mean - 1.96 * np.sqrt(f_var)
f_upper = f_mean + 1.96 * np.sqrt(f_var)
y_lower = y_mean - 1.96 * np.sqrt(y_var)
y_upper = y_mean + 1.96 * np.sqrt(y_var)


plt.plot(X, Y, "kx", mew=2, label="input data")
plt.plot(Xplot, f_mean, "-", color="C0", label="mean")
#plt.plot(Xplot, f_lower, "--", color="C0", label="f 95% confidence")
#plt.plot(Xplot, f_upper, "--", color="C0")
#plt.fill_between(
#    Xplot[:, 0], f_lower[:, 0], f_upper[:, 0], color="C0", alpha=0.1
#)
plt.plot(Xplot, y_lower, "--", color="C0", label="Y 95% confidence")
plt.plot(Xplot, y_upper, "--", color="C0")
plt.fill_between(
    Xplot[:, 0], y_lower[:, 0], y_upper[:, 0], color="C0", alpha=0.2
)
plt.legend()
plt.ylabel('Depth (z)')
plt.xlabel('Distance along Seafloor (x)')

gpflow.utilities.print_summary(model)

plt.show()