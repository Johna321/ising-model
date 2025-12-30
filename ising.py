import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from matplotlib.animation import FuncAnimation

class Ising:
    def __init__(self, L=8, J_meV=10):
        self.L = L
        self.J = J_meV
        self.k_B = 8.617333262e-2 # meV / K
        self.reset()
    
    def reset(self, T=400):
        self.rng = np.random.default_rng()
        self.T = T
        self.beta = 1.0 / (self.k_B * self.T)
        self.c = self.rng.choice([-1, 1], size=(self.L, self.L))
        self.energy_history = []
        self.mag_history = []
        self.sweep_count = 0
    
    # calculate the energy of an entire config
    def H(self):
        # kernel for 4 neighbors
        kernel = np.array([[0, 1, 0],
                            [1, 0, 1],
                            [0, 1, 0]])

        # sum of neighbors for each spin
        neighbor_sum = convolve2d(self.c, kernel, mode='same', boundary='fill', fillvalue=0)

        # spin-spin interactions on each, divide by 2 to avoid double count
        total_energy = -self.J * np.sum(self.c * neighbor_sum) / 2

        return total_energy

    # calculate energy delta between two configs (NO PERIODIC BOUNDARY)
    def local_H(self, i, j):
        # pad with zeros to avoid index errors
        padded = np.pad(self.c, pad_width=1, mode='constant', constant_values=0)

        # be mindful that our indices are -1 for the paddded config
        neighbors = [
            padded[i, j+1],   # up
            padded[i+2, j+1], # down
            padded[i+1, j],   # left
            padded[i+1, j+2]  # right
        ]

        return -self.J * self.c[i,j] * sum(neighbors)

    # metropolis-hastings D:
    def metropolis_sweep(self):
        for _ in range(self.L * self.L):
            # randomly select candidate spin for flip
            i = self.rng.integers(0,self.L)
            j = self.rng.integers(0,self.L)

            # adjust by factor of -2 to get true delta
            delta_H = -2 * self.local_H(i, j)

            # acceptance criteria
            if delta_H <= 0 or np.exp(-self.beta * delta_H) > self.rng.random():
                self.c[i,j] *= -1
            
        self.sweep_count += 1
        self.energy_history.append(self.H() / (self.L * self.L))
        self.mag_history.append(np.abs(np.mean(self.c)))

sim = Ising(L=150)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
ax_spin, ax_energy, ax_mag, ax_info = axes.flatten()

im = ax_spin.imshow(sim.c, cmap='RdBu', vmin=-1, vmax=1, interpolation='nearest')
ax_spin.set_title('configuration')
fig.colorbar(im, ax=ax_spin)

line_energy, = ax_energy.plot([], [], 'b-')
ax_energy.set_xlabel('sweep')
ax_energy.set_ylabel('E per spin')
ax_energy.set_title('energy')
ax_energy.grid(True)

line_mag, = ax_mag.plot([], [], 'r-')
ax_mag.set_xlabel('sweep')
ax_mag.set_ylabel('|mag|')
ax_mag.set_title('mag')
ax_mag.grid(True)

ax_energy.set_xlim(0, 150)
ax_energy.set_ylim(-30, 30)
ax_mag.set_xlim(0, 150)
ax_mag.set_ylim(0, 1.15)

ax_info.axis('off')
info_text = ax_info.text(0.1, 0.5, '')

def init():
    line_energy.set_data([], [])
    line_mag.set_data([], [])
    return [im, line_energy, line_mag, info_text]

def animate(frame):
    # run a sweep
    sim.metropolis_sweep()
    
    # update spin config
    im.set_data(sim.c)
    
    # update plots
    if len(sim.energy_history) > 0:
        sweeps = list(range(len(sim.energy_history)))
        line_energy.set_data(sweeps, sim.energy_history)
        line_mag.set_data(sweeps, sim.mag_history)
        
        
    
    # update info text
    info_text.set_text(
        f"sweep: {sim.sweep_count}\n"
        f"temp: {sim.T:.3f}\n"
        f"beta: {sim.beta:.3f}\n"
        f"energy: {sim.energy_history[-1] if sim.energy_history else 0:.3f}\n"
        f"mag: {sim.mag_history[-1] if sim.mag_history else 0:.3f}"
    )
    
    return [im, line_energy, line_mag, info_text]

# animation
anim = FuncAnimation(fig, animate, init_func=init, frames=100,
                    interval=50, blit=True, repeat=False)
                    
plt.tight_layout()
plt.show()

        
