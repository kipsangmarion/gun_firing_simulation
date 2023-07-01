# Bullet Firing Simulation

## Overview
This is a Python program that simulates the trajectory of a bullet fired from a gun at different angles. The simulation considers various physical factors including gravity, air resistance, wind speed and direction. The simulation also calculates the distance the bullet covered and the impact force when the bullet hits a stationary target.

## Variables and Constants

### Constants
1. `g`: Acceleration due to gravity (9.81 m/s²)
2. `air_density`: Air density at sea level and 15 degrees Celsius (1.225 kg/m³)

### Gun and Bullet Properties
1. `gun_height`: Height of the gun (1 m)
2. `initial_bullet_velocity`: Initial speed of the bullet (350 m/s)
3. `bullet_mass`: Mass of the bullet (0.007 kg)
4. `bullet_diameter`: Diameter of the bullet (0.009 m)
5. `bullet_drag_coefficient`: Drag coefficient of the bullet (0.05, dimensionless)

### Target Properties
1. `target_distance`: Horizontal distance to the target (500 m)
2. `target_size`: Diameter of the target (1 m)

### Wind Properties
1. `wind_speed`: Speed of the wind (10 m/s)
2. `wind_direction`: Direction of the wind (0 radians)

### Initial Conditions
1. `firing_angle`: Firing angle (45 degrees, converted to radians)
2. `x`: Initial horizontal position (0 m)
3. `y`: Initial vertical position (`gun_height` m)
4. `v`: Initial velocity (`initial_bullet_velocity` m/s)
5. `θ`: Initial angle (`firing_angle` radians)

### Time Settings
1. `t`: Initial time (0 s)
2. `dt`: Time step (0.01 s)

## Runge-Kutta Method
The Runge-Kutta method is used to numerically solve the differential equations of motion. 

## Results
The program outputs the trajectory of the bullet, the distance it covered, and if the bullet hit the target, the impact force. The impact force is calculated based on the bullet's momentum change during an assumed stopping time of 0.001 s.

## Running the Code
To run the code, simply run the script in a Python environment with `matplotlib` installed. You can adjust the variables and constants to see how they affect the bullet's trajectory and the resulting impact force. The plot will show the trajectory of the bullet.

## Notes
The simulation is a simplification and does not account for certain real-world factors like bullet spin, temperature variations, humidity, altitude, or changes in wind speed and direction over distance. The calculation of the impact force is also an overestimate and the actual force could be much less. More accurate simulation would require complex modeling and simulation, taking into account factors like deformation, energy dissipation, material properties, etc.
