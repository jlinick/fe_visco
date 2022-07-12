# Viscoelastic Finite Element Simulation of Supraglacial Lake Loading/Unloading

 
## Running the model

clone the respository

`git clone https://github.com/jlinick/fe_visco.git`

jump into the directory

`cd fe_visco`

run the simulation with
`./start_container.sh`

This command will build the docker image (if necessary), then jump you into a running container, call run_simulation.py to run the fe model, save products to a /products subdirectory, then exit and remove the docker container.

