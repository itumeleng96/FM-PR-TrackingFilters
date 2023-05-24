# Particle Filter

#### The Particle Filter consists of the following functions:
<li>createGaussianParticles() - creates random particles with each particle having a random initial state</li>
<li>predict() - using the process model and the previous state of the model, predict the next state of each particle</li>
<li>update() - compute updated weights of the particles based on their distance from the measured values </li>
<li>resample() - replace ineffictive particles with more effective particles</li>

