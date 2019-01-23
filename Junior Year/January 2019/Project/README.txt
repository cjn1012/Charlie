This OpenFoam folder simulates flow in a rectangular cross section as it instantly contracts into a smaller 'tube.'



To run this file, start blueCFD with the folder highlighted and run, in order:

1. blockMesh  - compiles files and generates a mesh that assists in plotting the final solution in OpenFoam
2. pimpleFoam - a large, timestep transient solver for an incompressible flow using PIMPLE algorithm, which is a PISO SIMPLE solver. 
3. paraFoam   - starts OpenFoam that allows you to view the results in a better user interface



Once paraFoam start OpenFoam up, there are a couple interesting filters that are interesting. To get to filters:

1. Press apply once OpenFoam starts
2. Press the play button (green triangle) at the top of the screen to load in the U, p, etc. data.
3. With the geometry showing on the screen, select U (velocity) at the top to see the velocity values in the steady state condition.

Now that you have the basic openFoam working, filters can be used to get helpful data and pictures detailing the flow. The first filter is called 'Glyph'.

Once you have selected glyph, it opens a new view within the pipeline browser. At first selection, the arrows detailing the direction of flow are huge. Scale them down
in the properties section using the bar.

